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
#   TOP_SEQS      if >0, fold only the TOP_SEQS best-scoring ProteinMPNN
#                 sequences per backbone (by MPNN score, lower = better);
#                 0 = fold all. Cuts production validation volume.
#   CAP           if >0, absolute cap on total YAMLs folded (keeps first CAP by
#                 sorted design id). Used by the smoke test; 0 = no cap.
set -eo pipefail

: "${MPNN_DIR:?set MPNN_DIR}"; : "${TARGET_SEQ:?set TARGET_SEQ}"
: "${OUT_DIR:?set OUT_DIR}"; : "${DIFF_SAMPLES:=1}"; : "${RECYCLES:=3}"
: "${TOP_SEQS:=0}"; : "${CAP:=0}"

mkdir -p "$OUT_DIR/yamls"

# Build one complex YAML per designed sequence (binder chain B + target chain A).
# Design id = <backbone>_s<si> where si is the sequence's ORIGINAL index in the
# ProteinMPNN fasta, so ids stay stable and join across Boltz/Chai/cross even
# when TOP_SEQS drops the lower-scoring sequences.
python - "$MPNN_DIR" "$OUT_DIR/yamls" "$TARGET_SEQ" "$TOP_SEQS" <<'PYEOF'
import glob, os, sys, re
mpnn_dir, yaml_dir, target_seq, top_seqs = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
n = 0
for fa in glob.glob(os.path.join(mpnn_dir, "seqs", "*.fa")):
    stem = os.path.basename(fa)[:-3]
    recs = open(fa).read().split(">")
    cands = []  # (score, si, seq)
    for si, rec in enumerate(recs[2:]):        # si = sample index -> unique per (backbone,sample)
        lines = rec.strip().split("\n")
        if len(lines) < 2:
            continue
        seq = "".join(lines[1:]).strip()
        if not seq:
            continue
        m = re.search(r"score=\s*([0-9.]+)", lines[0])   # ProteinMPNN score, lower is better
        score = float(m.group(1)) if m else 1e9
        cands.append((score, si, seq))
    if top_seqs > 0:
        cands = sorted(cands, key=lambda x: x[0])[:top_seqs]   # keep best-scoring per backbone
    for score, si, seq in cands:
        did = f"{stem}_s{si}".replace("/", "_").replace(".", "_")
        with open(os.path.join(yaml_dir, f"{did}.yaml"), "w") as o:
            o.write("version: 1\nsequences:\n")
            o.write(f"  - protein:\n      id: A\n      sequence: {target_seq}\n      msa: empty\n")
            o.write(f"  - protein:\n      id: B\n      sequence: {seq}\n      msa: empty\n")
        n += 1
print(f"[boltz] wrote {n} complex YAMLs (top_seqs={top_seqs})")
PYEOF

# Absolute cap (smoke test): keep only the first CAP YAMLs (SIGPIPE-safe).
if [ "$CAP" -gt 0 ]; then
  ls "$OUT_DIR"/yamls/*.yaml | sort | tail -n +$((CAP+1)) | xargs -r rm -f
  echo "[boltz] capped to $CAP complex YAMLs"
fi

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
