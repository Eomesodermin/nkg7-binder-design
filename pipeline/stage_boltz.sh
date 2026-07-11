#!/bin/bash
# stage_boltz.sh — Boltz-2 validation of designed binder:target complexes (one arm).
#
# For each ProteinMPNN sequence, build a binder:target complex YAML (single-seq)
# and co-fold with Boltz-2, emitting interface confidence (ipTM, complex pLDDT,
# per-chain pTM) that filter_designs.py consumes.
#
# We use Boltz-2 (pip env) NOT the ColabFold/Boltz-1 EasyBuild modules, which
# SIGILL on marvin's AVX-512-less AMD mlgpu nodes. Single-sequence mode
# (msa: empty) is standard for de novo binder validation and avoids the
# api.colabfold.com egress dependency on compute nodes.
#
# Required env vars:
#   MPNN_DIR      dir with seqs/*.fa (ProteinMPNN output)
#   TARGET_SEQ    target (receptor) amino-acid sequence
#   OUT_DIR       Boltz output dir
#   DIFF_SAMPLES  diffusion samples (1 smoketest, 5 production)
#   RECYCLES      recycling steps (3 default)
set -eo pipefail

: "${MPNN_DIR:?set MPNN_DIR}"; : "${TARGET_SEQ:?set TARGET_SEQ}"
: "${OUT_DIR:?set OUT_DIR}"; : "${DIFF_SAMPLES:=1}"; : "${RECYCLES:=3}"

mkdir -p "$OUT_DIR/yamls"

# Build one complex YAML per designed sequence (binder chain B + target chain A).
python - "$MPNN_DIR" "$OUT_DIR/yamls" "$TARGET_SEQ" <<'PYEOF'
import glob, os, sys
mpnn_dir, yaml_dir, target_seq = sys.argv[1], sys.argv[2], sys.argv[3]
n = 0
for fa in glob.glob(os.path.join(mpnn_dir, "seqs", "*.fa")):
    stem = os.path.basename(fa)[:-3]
    recs = open(fa).read().split(">")
    for si, rec in enumerate(recs[2:]):        # si = sample index -> unique per (backbone,sample)
        lines = rec.strip().split("\n")
        if len(lines) < 2:
            continue
        seq = "".join(lines[1:]).strip()
        if not seq:
            continue
        did = f"{stem}_s{si}".replace("/", "_").replace(".", "_")
        with open(os.path.join(yaml_dir, f"{did}.yaml"), "w") as o:
            o.write("version: 1\nsequences:\n")
            o.write(f"  - protein:\n      id: A\n      sequence: {target_seq}\n      msa: empty\n")
            o.write(f"  - protein:\n      id: B\n      sequence: {seq}\n      msa: empty\n")
        n += 1
print(f"[boltz] wrote {n} complex YAMLs")
PYEOF

echo "[boltz] predicting complexes (diffusion_samples=$DIFF_SAMPLES recycles=$RECYCLES)"
START=$SECONDS
# Boltz accepts a directory of YAMLs and folds each.
boltz predict "$OUT_DIR/yamls" \
  --out_dir "$OUT_DIR" \
  --recycling_steps "$RECYCLES" \
  --diffusion_samples "$DIFF_SAMPLES" \
  --no_kernels \
  --output_format pdb
echo "[boltz] runtime $((SECONDS-START))s"
NJ=$(find "$OUT_DIR" -name 'confidence_*_model_0.json' | wc -l)
echo "[boltz] scored complexes: $NJ"
[ "$NJ" -gt 0 ] && echo "STAGE_BOLTZ_OK"
