#!/usr/bin/env python
"""
render_fullscale.py — render SHARDED, conformer-aware SLURM job chains for the
full-scale multi-conformer NKG7 binder campaign.

Difference from render_jobs.py (smoke/production): here EVERY heavy stage is a
SLURM array (RFdiffusion, ProteinMPNN, AND Boltz — not just RFdiffusion), and
each arm is expanded into one independent chain PER CONFORMER (full-cross). A
chain directory is named "<arm>__<confid>" and holds rfdiffusion.sbatch,
proteinmpnn.sbatch, boltz.sbatch (all arrays), plus a "<arm>__consensus" dir
with chai.sbatch + crossspecies.sbatch (single, run once per arm on the pooled
top designs across all conformers).

Sharding contract (each array task computes its own input slice from
SLURM_ARRAY_TASK_ID, no pre-split step; outputs collision-free by design_id):
  RFd : 10000 bb / bb_per_shard(500) = 20 shards; shard k emits bb [k*500,(k+1)*500)
        via inference.design_startnum (distinct seed range + filenames).
  MPNN: sorted backbone list, shard k slices [k*500,(k+1)*500) backbones, parses
        just that slice, runs MPNN -> seqs/<bb>.fa (unique per backbone).
  Boltz: sorted backbone list, shard k slices [k*bb_per_boltz_shard,...) backbones,
         builds YAMLs for the top-N MPNN seqs of each, folds only its slice.
The wave-submission wrapper (submit_fullscale_waves.sh) submits these chains in
bounded waves so queued array tasks stay under the 300 MaxSubmit cap.
"""
import argparse, json, os, stat, yaml

ACCT = "ag_iei_abdullah"
CONDA = "source /opt/software/easybuild-INTEL/software/Miniforge3/24.1.2-0/etc/profile.d/conda.sh"

HDR = """#SBATCH --account={acct}
#SBATCH --partition={part}
#SBATCH --time={time}
{gpu}#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --job-name={jn}
{array}"""

def hdr(part, time, cpus, mem, jn, nshards=0, throttle=16, gpus=1):
    array = "#SBATCH --array=0-{0}%{1}".format(nshards-1, throttle) if nshards > 1 else ""
    gpu = "#SBATCH --gpus=1\n" if gpus else ""
    return HDR.format(acct=ACCT, part=part, time=time, cpus=cpus, mem=mem, jn=jn, array=array, gpu=gpu)


def rfd_sbatch(arm, conf, cfg, ws):
    d = cfg["defaults"]["rfdiffusion"]
    a = cfg["arms"][arm]
    n = d["num_designs_production"]
    per = d["bb_per_shard"]
    nshards = -(-n // per)              # ceil
    thr = d["array_throttle"]
    T = d["diffuser_T"]
    hs = "[" + ",".join(a["hotspot_res"]) + "]"
    rundir = f"{ws}/runs/fullscale/{arm}/{conf['id']}/01_rfdiffusion"
    cpdb = f"{ws}/{conf['pdb']}"
    tag = f"{arm}_{conf['id']}"
    return f"""#!/bin/bash
{hdr('mlgpu_short','08:00:00',8,'32G',f'rfd_{tag}',nshards,thr)}
set -eo pipefail
{CONDA}
conda activate SE3nv
module load CUDA/11.7.0
export DGLBACKEND=pytorch
mkdir -p {rundir}
SHARD=${{SLURM_ARRAY_TASK_ID:-0}}; NSHARDS={nshards}; PER={per}; N={n}
STARTNUM=$((SHARD*PER))
# last shard absorbs the remainder so the conformer gets exactly N backbones
if [ $SHARD -eq $((NSHARDS-1)) ]; then PER=$((N-STARTNUM)); fi
echo "[{tag}] RFd shard $SHARD/$NSHARDS: $PER designs from #$STARTNUM (target {conf['pdb']})"
START=$SECONDS
cd {ws}/software/RFdiffusion
python scripts/run_inference.py \\
  inference.output_prefix={rundir}/{tag} \\
  inference.input_pdb={cpdb} \\
  inference.num_designs=$PER \\
  inference.design_startnum=$STARTNUM \\
  inference.deterministic=True \\
  diffuser.T={T} \\
  'contigmap.contigs={a['contig']}' \\
  'ppi.hotspot_res={hs}'
echo "[{tag}] RFd shard $SHARD runtime $((SECONDS-START))s; bb now=$(ls {rundir}/{tag}_*.pdb 2>/dev/null | wc -l)"
echo STAGE_RFD_SHARD_OK
"""


def mpnn_sbatch(arm, conf, cfg, ws):
    d = cfg["defaults"]["proteinmpnn"]
    rd = cfg["defaults"]["rfdiffusion"]
    a = cfg["arms"][arm]
    n = rd["num_designs_production"]; per = d["bb_per_shard"]
    nshards = -(-n // per); thr = d["array_throttle"]
    nseq = d["num_seq_per_target_production"]
    sol = "--use_soluble_model" if d.get("use_soluble_model") else ""
    tchain = a["target_chain"]
    tag = f"{arm}_{conf['id']}"
    bbdir = f"{ws}/runs/fullscale/{arm}/{conf['id']}/01_rfdiffusion"
    outdir = f"{ws}/runs/fullscale/{arm}/{conf['id']}/02_proteinmpnn"
    return f"""#!/bin/bash
{hdr('mlgpu_short','04:00:00',8,'32G',f'mpnn_{tag}',nshards,thr)}
set -eo pipefail
{CONDA}
conda activate SE3nv
module load CUDA/11.7.0
MPNN={ws}/software/ProteinMPNN
SHARD=${{SLURM_ARRAY_TASK_ID:-0}}; PER={per}
SHARDDIR={outdir}/shard_$SHARD
mkdir -p {outdir}/seqs $SHARDDIR/in
# deterministic sorted backbone slice for this shard (avoid process substitution)
ls {bbdir}/{tag}_*.pdb | sort > $SHARDDIR/all_bb.txt
NBB=$(wc -l < $SHARDDIR/all_bb.txt)
START_I=$((SHARD*PER)); END_I=$((START_I+PER))
echo "[{tag}] MPNN shard $SHARD: backbones [$START_I,$END_I) of $NBB"
i=0
while read -r bb; do
  if [ $i -ge $START_I ] && [ $i -lt $END_I ]; then ln -sf "$bb" "$SHARDDIR/in/$(basename "$bb")"; fi
  i=$((i+1))
done < $SHARDDIR/all_bb.txt
if [ -z "$(ls -A $SHARDDIR/in 2>/dev/null)" ]; then echo "[{tag}] shard $SHARD empty, exit"; exit 0; fi
cd $MPNN
python helper_scripts/parse_multiple_chains.py --input_path=$SHARDDIR/in --output_path=$SHARDDIR/parsed.jsonl
FIRST=$(ls $SHARDDIR/in/*.pdb | head -1)
CHAINS=$(grep -h '^ATOM' $FIRST | cut -c22 | sort -u | tr -d ' ' | tr '\\n' ' ')
DESIGN=$(echo "$CHAINS" | tr ' ' '\\n' | grep -v "^{tchain}$" | tr '\\n' ' ')
python helper_scripts/assign_fixed_chains.py --input_path=$SHARDDIR/parsed.jsonl --output_path=$SHARDDIR/assigned.jsonl --chain_list "$DESIGN"
python protein_mpnn_run.py \\
  --jsonl_path $SHARDDIR/parsed.jsonl --chain_id_jsonl $SHARDDIR/assigned.jsonl \\
  --out_folder $SHARDDIR --num_seq_per_target {nseq} \\
  --sampling_temp "{d['sampling_temp']}" --model_name {d['model_name']} {sol} --batch_size 1
cp $SHARDDIR/seqs/*.fa {outdir}/seqs/ 2>/dev/null || true
echo "[{tag}] MPNN shard $SHARD produced $(ls $SHARDDIR/seqs/*.fa 2>/dev/null | wc -l) seq files"
echo STAGE_MPNN_SHARD_OK
"""


def boltz_sbatch(arm, conf, cfg, ws):
    d = cfg["defaults"]["validator"]
    rd = cfg["defaults"]["rfdiffusion"]
    a = cfg["arms"][arm]; sp = a["species"]
    mp = cfg["defaults"]["proteinmpnn"]
    ds = d["diffusion_samples_production"]
    top_seqs = mp["num_seq_per_target_production"]     # fold all N MPNN seqs (screen)
    n = rd["num_designs_production"]
    folds_per_shard = d["folds_per_shard"]
    bb_per_shard = max(1, folds_per_shard // top_seqs)  # backbones per boltz shard
    nshards = -(-n // bb_per_shard); thr = d["array_throttle"]
    tag = f"{arm}_{conf['id']}"
    mpnndir = f"{ws}/runs/fullscale/{arm}/{conf['id']}/02_proteinmpnn"
    outdir = f"{ws}/runs/fullscale/{arm}/{conf['id']}/03_boltz"
    tseq = f"{ws}/targets/{sp}_target.seq"
    return f"""#!/bin/bash
{hdr('mlgpu_short','08:00:00',8,'48G',f'boltz_{tag}',nshards,thr)}
set -eo pipefail
{CONDA}
conda activate boltz
SHARD=${{SLURM_ARRAY_TASK_ID:-0}}; BB_PER={bb_per_shard}
mkdir -p {outdir}/yamls
export MPNN_DIR="{mpnndir}"
export TARGET_SEQ="$(cat {tseq})"
export OUT_DIR="{outdir}"
export DIFF_SAMPLES={ds}
export TOP_SEQS={top_seqs}
export SHARD_IDX=$SHARD
export SHARD_BB={bb_per_shard}
# stage_boltz_sharded.sh: build YAMLs for backbones [SHARD*SHARD_BB,(SHARD+1)*SHARD_BB)
# (top TOP_SEQS MPNN seqs each) and fold only those with Boltz-2.
bash {ws}/pipeline/stage_boltz_sharded.sh
echo STAGE_BOLTZ_SHARD_OK
"""


def _refold_shards(arm, cfg, confids):
    """YAMLs/shard=250 (ds=5 -> 1250 folds/shard ~4h); array size covers the max shortlist."""
    BB = cfg["defaults"]["rfdiffusion"]["num_designs_production"]
    SEQ = cfg["defaults"]["proteinmpnn"]["num_seq_per_target_production"]
    frac = cfg["defaults"]["refold"]["shortlist_frac"]
    shortlist_max = int(len(confids) * BB * SEQ * frac)
    per = 250
    nshards = max(1, -(-shortlist_max // per))
    return per, nshards


def _consensus_dirs(arm, ws, confids):
    bl = [f"{ws}/runs/fullscale/{arm}/{c}/03_boltz" for c in confids]
    bb = [f"{ws}/runs/fullscale/{arm}/{c}/01_rfdiffusion" for c in confids]
    mp = [f"{ws}/runs/fullscale/{arm}/{c}/02_proteinmpnn" for c in confids]
    return bl, bb, mp


def _refold_common(arm, cfg, ws, confids):
    a = cfg["arms"][arm]; f = cfg["defaults"]["filters"]
    refdir = f"{ws}/runs/fullscale/{arm}/06_refold"
    return f"""export REFOLD_DIR="{refdir}"
export ARM="{arm}"
export FILTER_PY="{ws}/pipeline/filter_designs.py"
export BINDER_CHAIN="{a['binder_chain']}"; export TARGET_CHAIN="{a['target_chain']}"
export IPTM_MIN={f['iptm_min']}; export PLDDT_MIN={f['complex_plddt_min']}
export BPTM_MIN={f['binder_ptm_min']}; export RMSD_MAX={f['binder_rmsd_max']}"""


def refold_prep_sbatch(arm, cfg, ws, confids):
    """Single job: rank pooled ds=1 designs, shortlist top-frac, collect their YAMLs for ds=5."""
    frac = cfg["defaults"]["refold"]["shortlist_frac"]
    bl, bb, mp = _consensus_dirs(arm, ws, confids)
    return f"""#!/bin/bash
{hdr('intelsr_short','02:00:00',8,'32G',f'refprep_{arm}',gpus=0)}
set -eo pipefail
{CONDA}
conda activate boltz
{_refold_common(arm, cfg, ws, confids)}
export MODE=prep
export BOLTZ_DIR0="{bl[0]}"; export BB_DIR0="{bb[0]}"
export EXTRA_BOLTZ="{' '.join(bl[1:])}"; export EXTRA_BB="{' '.join(bb[1:])}"
export SHORTLIST_FRAC={frac}
bash {ws}/pipeline/stage_refold.sh
"""


def refold_fold_sbatch(arm, cfg, ws, confids):
    """Array: ds=5 re-fold the shortlisted YAMLs, 250 YAMLs/shard."""
    per, nshards = _refold_shards(arm, cfg, confids)
    thr = cfg["defaults"]["validator"]["array_throttle"]
    return f"""#!/bin/bash
{hdr('mlgpu_short','08:00:00',8,'48G',f'refold_{arm}',nshards,thr)}
set -eo pipefail
{CONDA}
conda activate boltz
{_refold_common(arm, cfg, ws, confids)}
export MODE=fold
export SHARD_IDX=${{SLURM_ARRAY_TASK_ID:-0}}
export SHARD_N={per}
bash {ws}/pipeline/stage_refold.sh
"""


def refold_rank_sbatch(arm, cfg, ws, confids):
    """Single job: re-rank the ds=5 folds on the MEAN of 5."""
    bl, bb, mp = _consensus_dirs(arm, ws, confids)
    return f"""#!/bin/bash
{hdr('intelsr_short','01:00:00',8,'32G',f'refrank_{arm}',gpus=0)}
set -eo pipefail
{CONDA}
conda activate boltz
{_refold_common(arm, cfg, ws, confids)}
export MODE=rank
export BB_DIR0="{bb[0]}"; export EXTRA_BB="{' '.join(bb[1:])}"
bash {ws}/pipeline/stage_refold.sh
echo STAGE_REFOLD_OK
"""


def chai_sbatch(arm, cfg, ws, confids):
    """Single (non-array) Chai consensus over the arm's ds=5 ranked shortlist (pooled across conformers)."""
    dcon = cfg["defaults"]["consensus"]; a = cfg["arms"][arm]; sp = a["species"]
    frac = dcon["topN_frac"]; tag = arm
    outdir = f"{ws}/runs/fullscale/{arm}/04_chai"
    ds5_ranked = f"{ws}/runs/fullscale/{arm}/06_refold/{arm}_ds5_ranked.csv"
    tseq = f"{ws}/targets/{sp}_target.seq"
    bl, bb, mp = _consensus_dirs(arm, ws, confids)
    mp0, mp_extra = mp[0], " ".join(mp[1:])
    return f"""#!/bin/bash
{hdr('mlgpu_short','08:00:00',8,'64G',f'chai_{tag}')}
set -eo pipefail
{CONDA}
conda activate chai
mkdir -p {outdir}
# Chai-fold the top {int(frac*100)}% of the ds=5-ranked shortlist (fold-stable, mean of 5).
# Absolute count computed at runtime from the ds=5 CSV row count (differs per arm).
NROWS=$(( $(wc -l < "{ds5_ranked}") - 1 ))
export CHAI_MAX=$(python3 -c "import math,sys; print(max(1, math.ceil($NROWS*{frac})))")
echo "[chai {arm}] ds5 shortlist=$NROWS -> top {int(frac*100)}% = $CHAI_MAX"
export CHAI_DOWNLOADS_DIR={ws}/software/chai_downloads
export MPNN_DIR="{mp0}"
export EXTRA_MPNN_DIRS="{mp_extra}"
export TARGET_SEQ="$(cat {tseq})"
export OUT_DIR="{outdir}"
export NATIVE_RANKED="{ds5_ranked}"
bash {ws}/pipeline/stage_chai.sh
echo STAGE_CHAI_OK
"""


def cross_sbatch(arm, cfg, ws, confids):
    dcr = cfg["defaults"]["crossspecies"]; a = cfg["arms"][arm]; sp = a["species"]
    cross_sp = cfg["species_map"][sp]; frac = dcr["topN_frac"]; tag = arm
    outdir = f"{ws}/runs/fullscale/{arm}/05_crossspecies"
    ds5_ranked = f"{ws}/runs/fullscale/{arm}/06_refold/{arm}_ds5_ranked.csv"
    cseq = f"{ws}/targets/{cross_sp}_target.seq"
    bl, bb, mp = _consensus_dirs(arm, ws, confids)
    mp0, mp_extra = mp[0], " ".join(mp[1:])
    return f"""#!/bin/bash
{hdr('mlgpu_short','08:00:00',8,'48G',f'cross_{tag}')}
set -eo pipefail
{CONDA}
conda activate boltz
mkdir -p {outdir}
# Cross-species re-fold the top {int(frac*100)}% of the ds=5-ranked shortlist vs the {cross_sp} target.
NROWS=$(( $(wc -l < "{ds5_ranked}") - 1 ))
export TOPN=$(python3 -c "import math,sys; print(max(1, math.ceil($NROWS*{frac})))")
echo "[cross {arm}] ds5 shortlist=$NROWS -> top {int(frac*100)}% = $TOPN"
export MPNN_DIR="{mp0}"
export EXTRA_MPNN_DIRS="{mp_extra}"
export NATIVE_RANKED="{ds5_ranked}"
export CROSS_SEQ="$(cat {cseq})"
export OUT_DIR="{outdir}"
export DIFF_SAMPLES={cfg['defaults']['refold']['diffusion_samples']}
bash {ws}/pipeline/stage_crossspecies.sh
echo STAGE_CROSS_OK
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--ws", required=True)
    ap.add_argument("--out", default="rendered_fullscale")
    ap.add_argument("--arms", nargs="*", default=None)
    args = ap.parse_args()
    cfg = yaml.safe_load(open(args.config))
    os.makedirs(args.out, exist_ok=True)
    arms = args.arms or list(cfg["arms"].keys())
    manifest = {"ws": args.ws, "chains": [], "arm_consensus": {}}
    for arm in arms:
        confs = cfg["arms"][arm]["conformers"]
        confids = [c["id"] for c in confs]
        for conf in confs:
            cdir = os.path.join(args.out, f"{arm}__{conf['id']}")
            os.makedirs(cdir, exist_ok=True)
            for stage, body in [("rfdiffusion", rfd_sbatch(arm, conf, cfg, args.ws)),
                                ("proteinmpnn", mpnn_sbatch(arm, conf, cfg, args.ws)),
                                ("boltz", boltz_sbatch(arm, conf, cfg, args.ws))]:
                p = os.path.join(cdir, f"{stage}.sbatch")
                open(p, "w").write(body); os.chmod(p, os.stat(p).st_mode | stat.S_IEXEC)
            manifest["chains"].append(f"{arm}__{conf['id']}")
        adir = os.path.join(args.out, f"{arm}__consensus")
        os.makedirs(adir, exist_ok=True)
        # ds=5 refold (prep -> fold-array -> rank) THEN chai + cross on the ds=5 shortlist.
        for stage, body in [("refold_prep", refold_prep_sbatch(arm, cfg, args.ws, confids)),
                            ("refold_fold", refold_fold_sbatch(arm, cfg, args.ws, confids)),
                            ("refold_rank", refold_rank_sbatch(arm, cfg, args.ws, confids)),
                            ("chai", chai_sbatch(arm, cfg, args.ws, confids)),
                            ("crossspecies", cross_sbatch(arm, cfg, args.ws, confids))]:
            p = os.path.join(adir, f"{stage}.sbatch")
            open(p, "w").write(body); os.chmod(p, os.stat(p).st_mode | stat.S_IEXEC)
        manifest["arm_consensus"][arm] = {"dir": f"{arm}__consensus",
            "stage_order": ["refold_prep", "refold_fold", "refold_rank", "chai", "crossspecies"],
            "depends_on_chains": [f"{arm}__{c}" for c in confids]}
    json.dump(manifest, open(os.path.join(args.out, "manifest_fullscale.json"), "w"), indent=2)
    print(f"rendered {len(manifest['chains'])} sharded chains + {len(arms)} consensus dirs -> {args.out}/")
    print(f"manifest: {args.out}/manifest_fullscale.json")


if __name__ == "__main__":
    main()
