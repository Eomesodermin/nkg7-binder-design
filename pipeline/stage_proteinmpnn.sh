#!/bin/bash
# stage_proteinmpnn.sh — inverse-fold RFdiffusion binder backbones (one arm).
#
# Designs sequences for the BINDER chain only, holding the target chain fixed.
# RFdiffusion writes the binder as the last chain; we detect it and fix all
# target-chain positions so only the binder is redesigned.
#
# Required env vars:
#   MPNN_DIR      path to ProteinMPNN clone
#   BACKBONE_DIR  dir of RFdiffusion *_N.pdb backbones
#   OUT_DIR       output dir for sequences
#   NUM_SEQ       sequences per backbone
#   SAMPLING_TEMP e.g. "0.1"
#   MODEL_NAME    v_48_020 (default; tolerant of RFdiffusion backbones)
#   USE_SOLUBLE   "1" to use the soluble-trained weights (expression-biased)
#   TARGET_CHAIN  target chain letter to hold FIXED (e.g. A); binder chain is redesigned
set -eo pipefail

: "${MPNN_DIR:?set MPNN_DIR}"; : "${BACKBONE_DIR:?set BACKBONE_DIR}"
: "${OUT_DIR:?set OUT_DIR}"; : "${NUM_SEQ:=8}"; : "${SAMPLING_TEMP:=0.1}"
: "${MODEL_NAME:=v_48_020}"; : "${USE_SOLUBLE:=1}"; : "${TARGET_CHAIN:=A}"

mkdir -p "$OUT_DIR"
cd "$MPNN_DIR"

SOL_FLAG=""; [ "$USE_SOLUBLE" = "1" ] && SOL_FLAG="--use_soluble_model"

# Parse chains from the first backbone: target chain = TARGET_CHAIN, the rest are binder.
FIRST=$(ls "$BACKBONE_DIR"/*.pdb | head -1)
CHAINS=$(grep '^ATOM' "$FIRST" | cut -c22 | sort -u | tr -d ' ' | tr '\n' ' ')
echo "[proteinmpnn] chains in backbone: $CHAINS ; fixing target chain $TARGET_CHAIN"

START=$SECONDS
# Build a parsed-chains + fixed-positions spec so ONLY the binder chain is designed.
python helper_scripts/parse_multiple_chains.py \
  --input_path="$BACKBONE_DIR" --output_path="$OUT_DIR/parsed_chains.jsonl"
# design chains = all chains except target; fixed chains = target
DESIGN_CHAINS=$(echo "$CHAINS" | tr ' ' '\n' | grep -v "^${TARGET_CHAIN}$" | tr '\n' ' ')
python helper_scripts/assign_fixed_chains.py \
  --input_path="$OUT_DIR/parsed_chains.jsonl" \
  --output_path="$OUT_DIR/assigned_chains.jsonl" \
  --chain_list "$DESIGN_CHAINS"

python protein_mpnn_run.py \
  --jsonl_path "$OUT_DIR/parsed_chains.jsonl" \
  --chain_id_jsonl "$OUT_DIR/assigned_chains.jsonl" \
  --out_folder "$OUT_DIR" \
  --num_seq_per_target "$NUM_SEQ" \
  --sampling_temp "$SAMPLING_TEMP" \
  --model_name "$MODEL_NAME" \
  $SOL_FLAG \
  --batch_size 1
echo "[proteinmpnn] runtime $((SECONDS-START))s"
N=$(ls "$OUT_DIR"/seqs/*.fa 2>/dev/null | wc -l)
echo "[proteinmpnn] produced sequence files: $N"
[ "$N" -gt 0 ] && echo "STAGE_MPNN_OK"
