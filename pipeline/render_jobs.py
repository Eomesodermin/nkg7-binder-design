#!/usr/bin/env python
"""
render_jobs.py — render per-arm SLURM job scripts for the binder pipeline from
a campaign config. Produces self-contained sbatch scripts (absolute paths
throughout) for each stage of each arm, plus a manifest the orchestrator uses
to submit and chain them.

This is the reusable core: point it at a different config (different target,
epitopes, hotspots) and it regenerates the whole campaign. Nothing here is
NKG7-specific.

Design notes baked into the rendered scripts (see DESIGN_DECISIONS.md):
  * Full-length target retained; hotspots restrict the interface to loops so the
    TM bundle acts as a steric wall.
  * RFdiffusion uses noise_scale=0 for sharper binder backbones.
  * ProteinMPNN fixes the target chain, redesigns only the binder, soluble model.
  * Boltz-2 single-sequence complex fold; filter on ipTM / complex pLDDT / binder ptm / design-vs-refold RMSD.

Usage:
  python render_jobs.py --config configs/nkg7_campaign.yaml \
      --mode smoketest --ws /lustre/.../nkg7_binder_design --out rendered/
"""
import argparse, os, yaml, json, stat

SBATCH_HEADER = """#SBATCH --account={account}
#SBATCH --partition={partition}
#SBATCH --time={time}
#SBATCH --gpus={gpus}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
"""

CONDA = "source /opt/software/easybuild-INTEL/software/Miniforge3/24.1.2-0/etc/profile.d/conda.sh"


def rfd_script(arm, cfg, ws, mode):
    d = cfg["defaults"]["rfdiffusion"]
    n = d["num_designs_smoketest"] if mode == "smoketest" else d["num_designs_production"]
    T = 25 if mode == "smoketest" else d["diffuser_T"]
    a = cfg["arms"][arm]
    hs = ",".join(a["hotspot_res"])
    rundir = f"{ws}/runs/{mode}/{arm}/01_rfdiffusion"
    target = f"{ws}/targets/{os.path.basename(a['target_pdb'])}"
    partition = "mlgpu_devel" if mode == "smoketest" else "mlgpu_short"
    time = "01:00:00" if mode == "smoketest" else "08:00:00"
    return f"""#!/bin/bash
{SBATCH_HEADER.format(account='ag_iei_abdullah', partition=partition, time=time, gpus=1, cpus=8, mem='32G')}
set -eo pipefail
{CONDA}
conda activate SE3nv
module load CUDA/11.7.0
export DGLBACKEND=pytorch
mkdir -p {rundir}
cd {ws}/software/RFdiffusion
echo "[{arm}] RFdiffusion N={n} T={T} hotspots={hs}"
START=$SECONDS
python scripts/run_inference.py \\
  inference.output_prefix={rundir}/{arm} \\
  inference.model_directory_path={ws}/software/RFdiffusion/models \\
  inference.input_pdb={target} \\
  inference.num_designs={n} \\
  'contigmap.contigs={a["contig"]}' \\
  'ppi.hotspot_res=[{hs}]' \\
  denoiser.noise_scale_ca={d["denoiser_noise_scale"]} \\
  denoiser.noise_scale_frame={d["denoiser_noise_scale"]} \\
  diffuser.T={T}
echo "[{arm}] RFdiffusion runtime $((SECONDS-START))s"
NB=$(ls {rundir}/{arm}_*.pdb 2>/dev/null | wc -l)
echo "[{arm}] produced $NB backbones"
[ "$NB" -gt 0 ] && echo "STAGE_RFD_OK"
"""


def mpnn_script(arm, cfg, ws, mode):
    d = cfg["defaults"]["proteinmpnn"]
    nseq = d["num_seq_per_target_smoketest"] if mode == "smoketest" else d["num_seq_per_target_production"]
    a = cfg["arms"][arm]
    bbdir = f"{ws}/runs/{mode}/{arm}/01_rfdiffusion"
    outdir = f"{ws}/runs/{mode}/{arm}/02_proteinmpnn"
    sol = "--use_soluble_model" if d["use_soluble_model"] else ""
    tchain = a["target_chain"]
    return f"""#!/bin/bash
{SBATCH_HEADER.format(account='ag_iei_abdullah', partition='mlgpu_devel' if mode=='smoketest' else 'mlgpu_short', time='00:30:00' if mode=='smoketest' else '02:00:00', gpus=1, cpus=8, mem='32G')}
set -eo pipefail
{CONDA}
conda activate SE3nv
module load CUDA/11.7.0
MPNN={ws}/software/ProteinMPNN
mkdir -p {outdir}
cd $MPNN
echo "[{arm}] ProteinMPNN nseq={nseq} (fix target chain {tchain}, redesign binder)"
START=$SECONDS
python helper_scripts/parse_multiple_chains.py --input_path={bbdir} --output_path={outdir}/parsed.jsonl
# design = all chains except the target chain
CHAINS=$(grep -h '^ATOM' {bbdir}/{arm}_0.pdb | cut -c22 | sort -u | tr -d ' ' | tr '\\n' ' ')
DESIGN=$(echo "$CHAINS" | tr ' ' '\\n' | grep -v "^{tchain}$" | tr '\\n' ' ')
echo "[{arm}] chains=$CHAINS design=$DESIGN"
python helper_scripts/assign_fixed_chains.py --input_path={outdir}/parsed.jsonl --output_path={outdir}/assigned.jsonl --chain_list "$DESIGN"
python protein_mpnn_run.py \\
  --jsonl_path {outdir}/parsed.jsonl \\
  --chain_id_jsonl {outdir}/assigned.jsonl \\
  --out_folder {outdir} \\
  --num_seq_per_target {nseq} \\
  --sampling_temp "{d['sampling_temp']}" \\
  --model_name {d['model_name']} {sol} --batch_size 1
echo "[{arm}] ProteinMPNN runtime $((SECONDS-START))s"
NS=$(ls {outdir}/seqs/*.fa 2>/dev/null | wc -l)
echo "[{arm}] produced $NS seq files"
[ "$NS" -gt 0 ] && echo "STAGE_MPNN_OK"
"""


def boltz_script(arm, cfg, ws, mode):
    d = cfg["defaults"]["validator"]
    ds = d["diffusion_samples_smoketest"] if mode == "smoketest" else d["diffusion_samples_production"]
    rec = d["recycling_steps"]
    cap = d.get("max_complexes_smoketest", 0) if mode == "smoketest" else d.get("max_complexes_production", 0)
    a = cfg["arms"][arm]
    sp = a["species"]
    mpnndir = f"{ws}/runs/{mode}/{arm}/02_proteinmpnn"
    outdir = f"{ws}/runs/{mode}/{arm}/03_boltz"
    yamldir = f"{outdir}/yamls"
    target_seq_file = f"{ws}/targets/{sp}_target.seq"
    return f"""#!/bin/bash
{SBATCH_HEADER.format(account='ag_iei_abdullah', partition='mlgpu_devel' if mode=='smoketest' else 'mlgpu_short', time='01:00:00' if mode=='smoketest' else '08:00:00', gpus=1, cpus=8, mem='48G')}
set -eo pipefail
{CONDA}
conda activate boltz
mkdir -p {yamldir} {outdir}
TARGET_SEQ=$(cat {target_seq_file})
# Build one binder:target complex YAML per ProteinMPNN design (single-seq mode).
python - "{mpnndir}" "{yamldir}" "$TARGET_SEQ" <<'PYEOF'
import glob, os, sys
mpnn_dir, yaml_dir, target_seq = sys.argv[1], sys.argv[2], sys.argv[3]
n=0
for fa in glob.glob(os.path.join(mpnn_dir,"seqs","*.fa")):
    stem=os.path.basename(fa)[:-3]
    recs=open(fa).read().split(">")
    for si,rec in enumerate(recs[2:]):   # si = sample index within this backbone -> unique
        lines=rec.strip().split("\\n")
        if len(lines)<2: continue
        seq="".join(lines[1:]).strip()
        if not seq: continue
        did=f"{{stem}}_s{{si}}".replace("/","_").replace(".","_")   # e.g. mouse_ecl2_7_s3
        with open(os.path.join(yaml_dir,f"{{did}}.yaml"),"w") as o:
            o.write("version: 1\\nsequences:\\n")
            o.write(f"  - protein:\\n      id: A\\n      sequence: {{target_seq}}\\n      msa: empty\\n")
            o.write(f"  - protein:\\n      id: B\\n      sequence: {{seq}}\\n      msa: empty\\n")
        n+=1
print(f"[boltz] wrote {{n}} complex YAMLs")
PYEOF
# Optional smoke-test cap: keep only the first N complex YAMLs.
CAP={cap}
if [ "$CAP" -gt 0 ]; then
  # SIGPIPE-safe cap: delete every YAML past the first $CAP (sorted).
  ls {yamldir}/*.yaml | sort | tail -n +$((CAP+1)) | xargs -r rm -f
  echo "[{arm}] capped to $CAP complexes for smoketest"
fi
echo "[{arm}] Boltz-2 folding complexes (diffusion_samples={ds} recycles={rec})"
START=$SECONDS
boltz predict {yamldir} --out_dir {outdir} --recycling_steps {rec} \\
  --diffusion_samples {ds} --no_kernels --output_format pdb
echo "[{arm}] Boltz runtime $((SECONDS-START))s"
NJ=$(find {outdir} -name 'confidence_*_model_0.json' | wc -l)
echo "[{arm}] scored $NJ complexes"
[ "$NJ" -gt 0 ] && echo "STAGE_BOLTZ_OK"
"""


def chai_script(arm, cfg, ws, mode):
    """Chai-1 orthogonal-consensus refold of the same binder:target complexes."""
    d = cfg["defaults"]["consensus"]
    v = cfg["defaults"]["validator"]
    nd = d["diffn_timesteps_smoketest"] if mode == "smoketest" else d["diffn_timesteps_production"]
    rec = d["recycling_steps"]
    cap = v.get("max_complexes_smoketest", 0) if mode == "smoketest" else v.get("max_complexes_production", 0)
    a = cfg["arms"][arm]
    sp = a["species"]
    mpnndir = f"{ws}/runs/{mode}/{arm}/02_proteinmpnn"
    outdir = f"{ws}/runs/{mode}/{arm}/04_chai"
    target_seq_file = f"{ws}/targets/{sp}_target.seq"
    partition = "mlgpu_devel" if mode == "smoketest" else "mlgpu_short"
    time = "01:00:00" if mode == "smoketest" else "08:00:00"
    return f"""#!/bin/bash
{SBATCH_HEADER.format(account='ag_iei_abdullah', partition=partition, time=time, gpus=1, cpus=8, mem='64G')}
set -eo pipefail
{CONDA}
conda activate chai
export CHAI_DOWNLOADS_DIR={ws}/software/chai_downloads
mkdir -p {outdir}
export MPNN_DIR="{mpnndir}"
export TARGET_SEQ="$(cat {target_seq_file})"
export OUT_DIR="{outdir}"
export CAP={cap}
export DIFFN_TIMESTEPS={nd}
export RECYCLES={rec}
bash {ws}/pipeline/stage_chai.sh
"""


def crossspecies_script(arm, cfg, ws, mode):
    """Fold this arm's top-N binders against the OTHER species' target (Boltz-2)."""
    d = cfg["defaults"]["crossspecies"]
    topn = d["topn_smoketest"] if mode == "smoketest" else d["topn_production"]
    ds = d["diffusion_samples"]
    rec = d["recycling_steps"]
    a = cfg["arms"][arm]
    sp = a["species"]
    cross_sp = d["species_map"][sp]
    mpnndir = f"{ws}/runs/{mode}/{arm}/02_proteinmpnn"
    outdir = f"{ws}/runs/{mode}/{arm}/05_crossspecies"
    native_ranked = f"{ws}/runs/{mode}/{arm}/04_filtered/{arm}_ranked.csv"
    cross_seq_file = f"{ws}/targets/{cross_sp}_target.seq"
    partition = "mlgpu_devel" if mode == "smoketest" else "mlgpu_short"
    time = "01:00:00" if mode == "smoketest" else "08:00:00"
    return f"""#!/bin/bash
{SBATCH_HEADER.format(account='ag_iei_abdullah', partition=partition, time=time, gpus=1, cpus=8, mem='48G')}
set -eo pipefail
{CONDA}
conda activate boltz
mkdir -p {outdir} {ws}/runs/{mode}/{arm}/04_filtered
# Rank the native-target designs first so we know the top-N to cross-fold.
python {ws}/pipeline/filter_designs.py \\
  --boltz-dir {ws}/runs/{mode}/{arm}/03_boltz \\
  --backbone-dir {ws}/runs/{mode}/{arm}/01_rfdiffusion \\
  --binder-chain {a['binder_chain']} --target-chain {a['target_chain']} \\
  --iptm-min {cfg['defaults']['filters']['iptm_min']} \\
  --plddt-min {cfg['defaults']['filters']['complex_plddt_min']} \\
  --binder-ptm-min {cfg['defaults']['filters']['binder_ptm_min']} \\
  --rmsd-max {cfg['defaults']['filters']['binder_rmsd_max']} \\
  --arm {arm} --out {native_ranked}
export MPNN_DIR="{mpnndir}"
export NATIVE_RANKED="{native_ranked}"
export CROSS_SEQ="$(cat {cross_seq_file})"
export OUT_DIR="{outdir}"
export TOPN={topn}
export DIFF_SAMPLES={ds}
export RECYCLES={rec}
bash {ws}/pipeline/stage_crossspecies.sh
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--mode", choices=["smoketest", "production"], default="smoketest")
    ap.add_argument("--ws", required=True, help="marvin workspace absolute path")
    ap.add_argument("--out", default="rendered")
    ap.add_argument("--arms", nargs="*", default=None, help="subset of arms (default: all)")
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config))
    os.makedirs(args.out, exist_ok=True)
    arms = args.arms or list(cfg["arms"].keys())
    manifest = {"mode": args.mode, "ws": args.ws, "arms": {}}
    for arm in arms:
        adir = os.path.join(args.out, args.mode, arm)
        os.makedirs(adir, exist_ok=True)
        scripts = {}
        stages = [("rfdiffusion", rfd_script), ("proteinmpnn", mpnn_script), ("boltz", boltz_script)]
        if cfg["defaults"].get("consensus", {}).get("enabled"):
            stages.append(("chai", chai_script))
        if cfg["defaults"].get("crossspecies", {}).get("enabled"):
            stages.append(("crossspecies", crossspecies_script))
        for stage, fn in stages:
            p = os.path.join(adir, f"{stage}.sbatch")
            with open(p, "w") as fh:
                fh.write(fn(arm, cfg, args.ws, args.mode))
            os.chmod(p, os.stat(p).st_mode | stat.S_IEXEC)
            scripts[stage] = p
        manifest["arms"][arm] = scripts
    with open(os.path.join(args.out, f"manifest_{args.mode}.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)
    nstages = len(next(iter(manifest["arms"].values()))) if manifest["arms"] else 0
    print(f"rendered {len(arms)} arms x {nstages} stages -> {args.out}/{args.mode}/")
    print(f"manifest: {args.out}/manifest_{args.mode}.json")


if __name__ == "__main__":
    main()
