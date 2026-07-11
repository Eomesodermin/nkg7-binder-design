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
    nshards = d.get("n_shards_smoketest", 1) if mode == "smoketest" else d.get("n_shards_production", 1)
    throttle = d.get("array_throttle", 2)
    a = cfg["arms"][arm]
    hs = ",".join(a["hotspot_res"])
    rundir = f"{ws}/runs/{mode}/{arm}/01_rfdiffusion"
    target = f"{ws}/targets/{os.path.basename(a['target_pdb'])}"
    partition = "mlgpu_devel" if mode == "smoketest" else "mlgpu_short"
    time = "01:00:00" if mode == "smoketest" else "08:00:00"
    # SLURM array: shard k runs designs [k*per, (k+1)*per) with a distinct seed,
    # numbering them contiguously via inference.design_startnum so files never collide.
    array_line = f"#SBATCH --array=0-{nshards-1}%{throttle}\n" if nshards > 1 else ""
    return f"""#!/bin/bash
{SBATCH_HEADER.format(account='ag_iei_abdullah', partition=partition, time=time, gpus=1, cpus=8, mem='32G')}{array_line}set -eo pipefail
{CONDA}
conda activate SE3nv
module load CUDA/11.7.0
export DGLBACKEND=pytorch
mkdir -p {rundir}
cd {ws}/software/RFdiffusion
SHARD=${{SLURM_ARRAY_TASK_ID:-0}}
NSHARDS={nshards}
PER=$(( {n} / NSHARDS ))
STARTNUM=$(( SHARD * PER ))
# last shard absorbs any remainder so the arm gets exactly {n} backbones
if [ "$SHARD" -eq "$(( NSHARDS - 1 ))" ]; then PER=$(( {n} - STARTNUM )); fi
# RFdiffusion has no inference.seed key; with inference.deterministic=True it
# seeds each design by its running index i_des. Distinct design_startnum per
# shard therefore gives each shard a distinct, non-overlapping seed range and
# unique output file numbers ({arm}_<startnum..>.pdb) — no collisions.
echo "[{arm}] RFdiffusion shard $SHARD/$NSHARDS: $PER designs from #$STARTNUM T={T} hotspots={hs}"
START=$SECONDS
python scripts/run_inference.py \\
  inference.output_prefix={rundir}/{arm} \\
  inference.model_directory_path={ws}/software/RFdiffusion/models \\
  inference.input_pdb={target} \\
  inference.num_designs=$PER \\
  inference.design_startnum=$STARTNUM \\
  inference.deterministic=True \\
  'contigmap.contigs={a["contig"]}' \\
  'ppi.hotspot_res=[{hs}]' \\
  denoiser.noise_scale_ca={d["denoiser_noise_scale"]} \\
  denoiser.noise_scale_frame={d["denoiser_noise_scale"]} \\
  diffuser.T={T}
echo "[{arm}] RFdiffusion shard $SHARD runtime $((SECONDS-START))s"
NB=$(ls {rundir}/{arm}_*.pdb 2>/dev/null | wc -l)
echo "[{arm}] backbones present after shard $SHARD: $NB"
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
    # production: fold only the top-N ProteinMPNN sequences per backbone (by score)
    top_seqs = d.get("top_seqs_per_backbone_smoketest", 0) if mode == "smoketest" else d.get("top_seqs_per_backbone_production", 0)
    a = cfg["arms"][arm]
    sp = a["species"]
    mpnndir = f"{ws}/runs/{mode}/{arm}/02_proteinmpnn"
    outdir = f"{ws}/runs/{mode}/{arm}/03_boltz"
    yamldir = f"{outdir}/yamls"
    target_seq_file = f"{ws}/targets/{sp}_target.seq"
    return f"""#!/bin/bash
{SBATCH_HEADER.format(account='ag_iei_abdullah', partition='mlgpu_devel' if mode=='smoketest' else 'mlgpu_short', time='01:00:00' if mode=='smoketest' else '08:00:00', gpus=1, cpus=8, mem='48G')}set -eo pipefail
{CONDA}
conda activate boltz
mkdir -p {yamldir} {outdir}
export MPNN_DIR="{mpnndir}"
export TARGET_SEQ="$(cat {target_seq_file})"
export OUT_DIR="{outdir}"
export DIFF_SAMPLES={ds}
export RECYCLES={rec}
export TOP_SEQS={top_seqs}
export CAP={cap}
# stage_boltz.sh builds YAMLs (top-N seqs/backbone if TOP_SEQS>0, then absolute
# CAP if >0) and folds them with Boltz-2.
bash {ws}/pipeline/stage_boltz.sh
"""


def chai_script(arm, cfg, ws, mode):
    """Chai-1 orthogonal-consensus refold of the TOP Boltz designs."""
    d = cfg["defaults"]["consensus"]
    v = cfg["defaults"]["validator"]
    f = cfg["defaults"]["filters"]
    nd = d["diffn_timesteps_smoketest"] if mode == "smoketest" else d["diffn_timesteps_production"]
    rec = d["recycling_steps"]
    # smoke: mirror the small Boltz cap; production: fold the top-N Boltz designs.
    chai_max = v.get("max_complexes_smoketest", 0) if mode == "smoketest" else v.get("chai_max_complexes_production", 0)
    a = cfg["arms"][arm]
    sp = a["species"]
    mpnndir = f"{ws}/runs/{mode}/{arm}/02_proteinmpnn"
    outdir = f"{ws}/runs/{mode}/{arm}/04_chai"
    ranked = f"{outdir}/{arm}_boltz_ranked.csv"   # own copy — avoids racing crossspecies
    target_seq_file = f"{ws}/targets/{sp}_target.seq"
    partition = "mlgpu_devel" if mode == "smoketest" else "mlgpu_short"
    time = "01:00:00" if mode == "smoketest" else "08:00:00"
    return f"""#!/bin/bash
{SBATCH_HEADER.format(account='ag_iei_abdullah', partition=partition, time=time, gpus=1, cpus=8, mem='64G')}set -eo pipefail
{CONDA}
conda activate boltz
mkdir -p {outdir}
# Rank the native Boltz designs so Chai folds the top ones (consensus on the shortlist).
python {ws}/pipeline/filter_designs.py \\
  --boltz-dir {ws}/runs/{mode}/{arm}/03_boltz \\
  --backbone-dir {ws}/runs/{mode}/{arm}/01_rfdiffusion \\
  --binder-chain {a['binder_chain']} --target-chain {a['target_chain']} \\
  --iptm-min {f['iptm_min']} --plddt-min {f['complex_plddt_min']} \\
  --binder-ptm-min {f['binder_ptm_min']} --rmsd-max {f['binder_rmsd_max']} \\
  --arm {arm} --out {ranked}
conda activate chai
export CHAI_DOWNLOADS_DIR={ws}/software/chai_downloads
export MPNN_DIR="{mpnndir}"
export TARGET_SEQ="$(cat {target_seq_file})"
export OUT_DIR="{outdir}"
export NATIVE_RANKED="{ranked}"
export CHAI_MAX={chai_max}
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
    native_ranked = f"{outdir}/{arm}_boltz_ranked.csv"   # own copy — no cross-stage race
    cross_seq_file = f"{ws}/targets/{cross_sp}_target.seq"
    partition = "mlgpu_devel" if mode == "smoketest" else "mlgpu_short"
    time = "01:00:00" if mode == "smoketest" else "08:00:00"
    return f"""#!/bin/bash
{SBATCH_HEADER.format(account='ag_iei_abdullah', partition=partition, time=time, gpus=1, cpus=8, mem='48G')}set -eo pipefail
{CONDA}
conda activate boltz
mkdir -p {outdir}
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
