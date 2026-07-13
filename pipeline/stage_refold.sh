#!/bin/bash
# stage_refold.sh — ds=5 wide-shortlist re-fold for the full-scale campaign.
#
# The ds=1 screen ranks designs cheaply but noisily (Spearman prod-vs-ds5 ~0.47).
# This stage re-folds the top SHORTLIST_FRAC of the pooled ds=1 designs at ds=5
# and re-ranks on the MEAN of 5 folds, so Chai/cross downstream consume a
# fold-stable shortlist rather than the noisy ds=1 point estimate.
#
# Three modes (MODE env):
#   prep : pool + rank ds=1 designs across all conformers (filter_designs.py),
#          take the top SHORTLIST_FRAC, and COLLECT their already-built ds=1
#          YAMLs into REFOLD_DIR/yamls/ (re-folding the same complex at ds=5 —
#          no sequence rebuild needed).
#   fold : array task; boltz-predict REFOLD_DIR/yamls slice [SHARD*SHARD_N,...)
#          at DIFF_SAMPLES=5 into REFOLD_DIR/shard_<SHARD>.
#   rank : filter_designs.py --aggregate mean over REFOLD_DIR -> ds5 ranked CSV.
#
# Required env (all modes): REFOLD_DIR, ARM, BINDER_CHAIN, TARGET_CHAIN,
#   IPTM_MIN, PLDDT_MIN, BPTM_MIN, RMSD_MAX, FILTER_PY
# prep: BOLTZ_DIR0 BB_DIR0 EXTRA_BOLTZ EXTRA_BB SHORTLIST_FRAC
# fold: SHARD_IDX SHARD_N   (+ boltz env)
# rank: STABILITY_PENALTY (optional)
set -eo pipefail
: "${MODE:?set MODE=prep|fold|rank}"; : "${REFOLD_DIR:?}"; : "${ARM:?}"; : "${FILTER_PY:?}"
: "${BINDER_CHAIN:=B}"; : "${TARGET_CHAIN:=A}"
: "${IPTM_MIN:=0.5}"; : "${PLDDT_MIN:=0.7}"; : "${BPTM_MIN:=0.5}"; : "${RMSD_MAX:=2.0}"

DS1_RANKED="$REFOLD_DIR/ds1_pooled_ranked.csv"
DS5_RANKED="$REFOLD_DIR/${ARM}_ds5_ranked.csv"
mkdir -p "$REFOLD_DIR/yamls"

if [ "$MODE" = "prep" ]; then
  : "${BOLTZ_DIR0:?}"; : "${BB_DIR0:?}"; : "${SHORTLIST_FRAC:=0.05}"
  EXTRA_BOLTZ_ARG=""; [ -n "${EXTRA_BOLTZ:-}" ] && EXTRA_BOLTZ_ARG="--extra-boltz-dirs ${EXTRA_BOLTZ}"
  EXTRA_BB_ARG="";    [ -n "${EXTRA_BB:-}" ]    && EXTRA_BB_ARG="--extra-backbone-dirs ${EXTRA_BB}"
  echo "[refold-prep $ARM] ranking ds=1 designs pooled across conformers"
  python "$FILTER_PY" \
    --boltz-dir "$BOLTZ_DIR0" --backbone-dir "$BB_DIR0" \
    --binder-chain "$BINDER_CHAIN" --target-chain "$TARGET_CHAIN" \
    --iptm-min "$IPTM_MIN" --plddt-min "$PLDDT_MIN" \
    --binder-ptm-min "$BPTM_MIN" --rmsd-max "$RMSD_MAX" \
    --arm "$ARM" --out "$DS1_RANKED" $EXTRA_BOLTZ_ARG $EXTRA_BB_ARG
  # take top SHORTLIST_FRAC by rank (CSV is already composite-ranked; row 1 = header)
  python - "$DS1_RANKED" "$SHORTLIST_FRAC" "$REFOLD_DIR/shortlist_ids.txt" <<'PYEOF'
import sys, csv, math
ranked, frac, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]
rows = list(csv.DictReader(open(ranked)))
k = max(1, math.ceil(len(rows) * frac))
ids = [r["design"] for r in rows[:k]]
open(out, "w").write("\n".join(ids) + "\n")
print(f"[refold-prep] {len(rows)} ds=1 designs -> top {k} ({frac*100:.0f}%) shortlisted")
PYEOF
  # collect each shortlisted design's existing ds=1 YAML (from any conformer's shard dirs)
  n=0
  while read -r did; do
    [ -z "$did" ] && continue
    y=$(find $BOLTZ_DIR0 ${EXTRA_BOLTZ:-} -type f -name "${did}.yaml" 2>/dev/null | head -1)
    if [ -n "$y" ]; then cp "$y" "$REFOLD_DIR/yamls/"; n=$((n+1)); fi
  done < "$REFOLD_DIR/shortlist_ids.txt"
  echo "[refold-prep $ARM] collected $n YAMLs for ds=5 re-fold -> $REFOLD_DIR/yamls"
  echo STAGE_REFOLD_PREP_OK

elif [ "$MODE" = "fold" ]; then
  : "${SHARD_IDX:?}"; : "${SHARD_N:?}"
  SHARD_OUT="$REFOLD_DIR/shard_${SHARD_IDX}"; mkdir -p "$SHARD_OUT/yamls"
  ls "$REFOLD_DIR"/yamls/*.yaml 2>/dev/null | sort > "$SHARD_OUT/all.txt" || true
  NY=$(wc -l < "$SHARD_OUT/all.txt")
  S=$((SHARD_IDX*SHARD_N)); E=$((S+SHARD_N))
  echo "[refold-fold $ARM] shard $SHARD_IDX: YAMLs [$S,$E) of $NY at ds=5"
  if [ "$S" -ge "$NY" ]; then echo "[refold-fold $ARM] empty shard"; echo STAGE_REFOLD_FOLD_OK; exit 0; fi
  i=0; while read -r y; do
    if [ $i -ge $S ] && [ $i -lt $E ]; then cp "$y" "$SHARD_OUT/yamls/"; fi
    i=$((i+1)); done < "$SHARD_OUT/all.txt"
  boltz predict "$SHARD_OUT/yamls" --out_dir "$SHARD_OUT" \
    --recycling_steps "${RECYCLES:=3}" --diffusion_samples 5 --no_kernels --output_format pdb
  echo "[refold-fold $ARM] shard $SHARD_IDX scored $(find $SHARD_OUT -name 'confidence_*_model_0.json'|wc -l)"
  echo STAGE_REFOLD_FOLD_OK

elif [ "$MODE" = "rank" ]; then
  echo "[refold-rank $ARM] re-ranking ds=5 folds on MEAN of 5"
  python "$FILTER_PY" \
    --boltz-dir "$REFOLD_DIR" --backbone-dir "$BB_DIR0" \
    --binder-chain "$BINDER_CHAIN" --target-chain "$TARGET_CHAIN" \
    --iptm-min "$IPTM_MIN" --plddt-min "$PLDDT_MIN" \
    --binder-ptm-min "$BPTM_MIN" --rmsd-max "$RMSD_MAX" \
    --arm "$ARM" --aggregate mean --out "$DS5_RANKED" ${EXTRA_BB:+--extra-backbone-dirs $EXTRA_BB}
  echo "[refold-rank $ARM] ds=5 ranked -> $DS5_RANKED ($(tail -n +2 $DS5_RANKED|wc -l) designs)"
  echo STAGE_REFOLD_RANK_OK
else
  echo "unknown MODE=$MODE"; exit 1
fi
