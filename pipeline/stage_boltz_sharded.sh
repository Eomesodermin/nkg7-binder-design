#!/bin/bash
# stage_boltz_sharded.sh — SHARDED Boltz-2 validation for the full-scale campaign.
#
# Same folding contract as stage_boltz.sh, but each SLURM array task processes
# only its slice of backbones: sorted seq-file list, this shard takes indices
# [SHARD_IDX*SHARD_BB, (SHARD_IDX+1)*SHARD_BB). YAMLs + Boltz outputs go to a
# shard-specific subdir so parallel shards never collide; a downstream harvest
# pools all shard_* subdirs.
#
# Required env vars (set by render_fullscale.py boltz_sbatch):
#   MPNN_DIR   dir with seqs/*.fa (pooled across this conformer's MPNN shards)
#   TARGET_SEQ target amino-acid sequence
#   OUT_DIR    boltz output dir for this conformer (03_boltz)
#   DIFF_SAMPLES  diffusion samples (1 for the ds=1 screen)
#   TOP_SEQS   fold the top-N MPNN seqs/backbone (5 for the screen)
#   SHARD_IDX  this array task's index
#   SHARD_BB   backbones per shard
set -eo pipefail
: "${MPNN_DIR:?set MPNN_DIR}"; : "${TARGET_SEQ:?set TARGET_SEQ}"; : "${OUT_DIR:?set OUT_DIR}"
: "${DIFF_SAMPLES:=1}"; : "${RECYCLES:=3}"; : "${TOP_SEQS:=5}"
: "${SHARD_IDX:?set SHARD_IDX}"; : "${SHARD_BB:?set SHARD_BB}"

SHARD_OUT="$OUT_DIR/shard_${SHARD_IDX}"
mkdir -p "$SHARD_OUT/yamls"

# Deterministic sorted backbone (seq-file) slice for this shard.
ls "$MPNN_DIR"/seqs/*.fa 2>/dev/null | sort > "$SHARD_OUT/all_seqs.txt" || true
NSEQFILES=$(wc -l < "$SHARD_OUT/all_seqs.txt")
START_I=$((SHARD_IDX*SHARD_BB)); END_I=$((START_I+SHARD_BB))
echo "[boltz-shard $SHARD_IDX] seq-files [$START_I,$END_I) of $NSEQFILES"
if [ "$START_I" -ge "$NSEQFILES" ]; then echo "[boltz-shard $SHARD_IDX] no backbones, exit"; echo STAGE_BOLTZ_SHARD_OK; exit 0; fi

# Build YAMLs for this shard's backbones only (top-TOP_SEQS MPNN seqs each).
python - "$SHARD_OUT/all_seqs.txt" "$SHARD_OUT/yamls" "$TARGET_SEQ" "$TOP_SEQS" "$START_I" "$END_I" <<'PYEOF'
import os, sys, re
listfile, yaml_dir, target_seq, top_seqs, start_i, end_i = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
fas = [l.strip() for l in open(listfile) if l.strip()][start_i:end_i]
n = 0
for fa in fas:
    stem = os.path.basename(fa)[:-3]
    recs = open(fa).read().split(">")
    cands = []
    for si, rec in enumerate(recs[2:]):          # skip fasta header + native seq
        lines = rec.strip().split("\n")
        if len(lines) < 2: continue
        seq = "".join(lines[1:]).strip()
        if not seq: continue
        m = re.search(r"score=\s*([0-9.]+)", lines[0])
        cands.append((float(m.group(1)) if m else 1e9, si, seq))
    if top_seqs > 0:
        cands = sorted(cands, key=lambda x: x[0])[:top_seqs]
    for score, si, seq in cands:
        did = f"{stem}_s{si}".replace("/", "_").replace(".", "_")
        with open(os.path.join(yaml_dir, f"{did}.yaml"), "w") as o:
            o.write("version: 1\nsequences:\n")
            o.write(f"  - protein:\n      id: A\n      sequence: {target_seq}\n      msa: empty\n")
            o.write(f"  - protein:\n      id: B\n      sequence: {seq}\n      msa: empty\n")
        n += 1
print(f"[boltz-shard] wrote {n} complex YAMLs from {len(fas)} backbones (top_seqs={top_seqs})")
PYEOF

echo "[boltz-shard $SHARD_IDX] folding (diffusion_samples=$DIFF_SAMPLES)"
START=$SECONDS
boltz predict "$SHARD_OUT/yamls" \
  --out_dir "$SHARD_OUT" \
  --recycling_steps "$RECYCLES" \
  --diffusion_samples "$DIFF_SAMPLES" \
  --no_kernels \
  --output_format pdb
echo "[boltz-shard $SHARD_IDX] runtime $((SECONDS-START))s"
NJ=$(find "$SHARD_OUT" -name 'confidence_*_model_0.json' | wc -l)
echo "[boltz-shard $SHARD_IDX] scored complexes: $NJ"
[ "$NJ" -gt 0 ] && echo STAGE_BOLTZ_SHARD_OK
