#!/bin/bash
# stage_refine.sh — Round-2 partial-diffusion refinement (one arm).
#
# Implements "Option I" from Tan et al. 2026 (bioRxiv 2026.03.04.709551):
# take a winning binder backbone, partially noise + denoise it with RFdiffusion
# (diffuser.partial_T) to generate a diverse family of related backbones around
# the same docked pose, then hand them to ProteinMPNN for re-sequencing. Used to
# lift fold quality (complex pLDDT) on interface-confident (high-ipTM) hits
# without discarding the binding mode.
#
# The target (chain A, 165 res) is held FIXED as a motif (contig "A1-165"); only
# the binder (chain B, length L) is resampled. partial_T controls how far the
# backbone is noised: small = gentle polish, large = more diversity / more drift.
#
# Required env vars:
#   RFD_DIR        path to RFdiffusion install
#   MODELS_DIR     path to RFdiffusion model weights
#   SEED_DIR       dir holding the round-1 seed backbone PDBs (target A + binder B)
#   SEED_LIST      space-separated backbone stems, e.g. "mouse_ecl1_355 mouse_ecl1_2"
#   OUT_PREFIX     output path prefix for refined backbones
#   PARTIAL_T      partial diffusion timesteps (e.g. 5 / 10 / 20)
#   VARIANTS       number of refined backbones per seed
#   TARGET_LEN     target chain length (default 165 for NKG7)
#   DIFFUSER_T     full trajectory length the model was trained on (default 50)
# Optional:
#   NOISE_SCALE    denoiser noise (default 0 -> sharper backbones)
set -eo pipefail

: "${RFD_DIR:?set RFD_DIR}"; : "${MODELS_DIR:?set MODELS_DIR}"
: "${SEED_DIR:?set SEED_DIR}"; : "${SEED_LIST:?set SEED_LIST}"
: "${OUT_PREFIX:?set OUT_PREFIX}"; : "${PARTIAL_T:?set PARTIAL_T}"
: "${VARIANTS:=20}"; : "${TARGET_LEN:=165}"; : "${DIFFUSER_T:=50}"
: "${NOISE_SCALE:=0}"

mkdir -p "$(dirname "$OUT_PREFIX")"
cd "$RFD_DIR"

total=0
for seed in $SEED_LIST; do
  seed_pdb="$SEED_DIR/${seed}.pdb"
  if [ ! -f "$seed_pdb" ]; then echo "[refine] MISSING seed $seed_pdb — skip"; continue; fi
  # binder length = residues in chain B of the seed PDB
  L=$(grep -E '^ATOM' "$seed_pdb" | awk '$5=="B"{print $6}' | sort -nu | wc -l)
  if [ "$L" -lt 1 ]; then echo "[refine] no chain B in $seed_pdb — skip"; continue; fi
  contig="[A1-${TARGET_LEN}/0 ${L}-${L}]"
  echo "[refine] seed=$seed L=$L contig=$contig partial_T=$PARTIAL_T variants=$VARIANTS"
  START=$SECONDS
  python scripts/run_inference.py \
    inference.output_prefix="${OUT_PREFIX}_${seed}_pt${PARTIAL_T}" \
    inference.model_directory_path="$MODELS_DIR" \
    inference.input_pdb="$seed_pdb" \
    inference.num_designs="$VARIANTS" \
    inference.deterministic=True \
    "contigmap.contigs=$contig" \
    denoiser.noise_scale_ca="$NOISE_SCALE" \
    denoiser.noise_scale_frame="$NOISE_SCALE" \
    diffuser.T="$DIFFUSER_T" \
    diffuser.partial_T="$PARTIAL_T"
  n=$(ls "${OUT_PREFIX}_${seed}_pt${PARTIAL_T}"_*.pdb 2>/dev/null | wc -l)
  echo "[refine] seed=$seed produced $n backbones ($((SECONDS-START))s)"
  total=$((total+n))
done
echo "[refine] TOTAL refined backbones: $total"
[ "$total" -gt 0 ] && echo "STAGE_REFINE_OK"
