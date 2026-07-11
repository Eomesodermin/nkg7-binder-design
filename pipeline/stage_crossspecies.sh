#!/bin/bash
# stage_crossspecies.sh — cross-species reactivity check (one arm).
#
# Folds this arm's top-N designed binders against the OTHER species' target
# with Boltz-2 (the same validator, so scores are directly comparable to the
# native-target run). A binder that keeps ipTM > threshold against BOTH its
# native and the cross target is predicted cross-reactive — valuable for a
# single reagent that works in human and mouse.
#
# Design ids are preserved (<stem>_s<si>) so cross scores join the native
# ranked table on `design`.
#
# Required env vars:
#   MPNN_DIR       dir with seqs/*.fa (native arm's ProteinMPNN output)
#   NATIVE_RANKED  native-target ranked CSV for this arm (picks top-N by score)
#   CROSS_SEQ      the OTHER species' target amino-acid sequence
#   OUT_DIR        cross-species Boltz output dir (e.g. .../05_crossspecies)
#   TOPN           number of top designs to cross-fold (default 10)
#   DIFF_SAMPLES   diffusion samples (default 1); RECYCLES (default 3)
set -eo pipefail

: "${MPNN_DIR:?set MPNN_DIR}"; : "${NATIVE_RANKED:?set NATIVE_RANKED}"
: "${CROSS_SEQ:?set CROSS_SEQ}"; : "${OUT_DIR:?set OUT_DIR}"
: "${TOPN:=10}"; : "${DIFF_SAMPLES:=1}"; : "${RECYCLES:=3}"

mkdir -p "$OUT_DIR/yamls"

# Build cross-target complex YAMLs for the top-N designs. We recover each
# design's binder sequence from the ProteinMPNN seqs using its design id.
python - "$MPNN_DIR" "$NATIVE_RANKED" "$CROSS_SEQ" "$OUT_DIR/yamls" "$TOPN" <<'PYEOF'
import glob, os, sys, csv
mpnn_dir, ranked, cross_seq, yaml_dir, topn = sys.argv[1:6]
topn = int(topn)

# map design id -> binder sequence (same parsing as stage_boltz)
seqmap = {}
for fa in glob.glob(os.path.join(mpnn_dir, "seqs", "*.fa")):
    stem = os.path.basename(fa)[:-3]
    recs = open(fa).read().split(">")
    for si, rec in enumerate(recs[2:]):
        lines = rec.strip().split("\n")
        if len(lines) < 2: continue
        seq = "".join(lines[1:]).strip()
        if not seq: continue
        did = f"{stem}_s{si}".replace("/", "_").replace(".", "_")
        seqmap[did] = seq

rows = list(csv.DictReader(open(ranked)))
# native ranked is already sorted best-first; take top-N with a known sequence
picked, n = [], 0
for r in rows:
    d = r["design"]
    if d in seqmap:
        picked.append(d); n += 1
    if n >= topn: break

for did in picked:
    with open(os.path.join(yaml_dir, f"{did}.yaml"), "w") as o:
        o.write("version: 1\nsequences:\n")
        o.write(f"  - protein:\n      id: A\n      sequence: {cross_seq}\n      msa: empty\n")
        o.write(f"  - protein:\n      id: B\n      sequence: {seqmap[did]}\n      msa: empty\n")
print(f"[cross] wrote {len(picked)} cross-target YAMLs (top {topn})")
PYEOF

echo "[cross] folding top designs against cross-species target (Boltz-2)"
START=$SECONDS
boltz predict "$OUT_DIR/yamls" \
  --out_dir "$OUT_DIR" \
  --recycling_steps "$RECYCLES" \
  --diffusion_samples "$DIFF_SAMPLES" \
  --no_kernels \
  --output_format pdb
echo "[cross] runtime $((SECONDS-START))s"
NJ=$(find "$OUT_DIR" -name 'confidence_*_model_0.json' | wc -l)
echo "[cross] scored cross complexes: $NJ"
[ "$NJ" -gt 0 ] && echo "STAGE_CROSS_OK"
