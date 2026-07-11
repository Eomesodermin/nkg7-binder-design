#!/bin/bash
# stage_rfdiffusion.sh — RFdiffusion binder-backbone generation (one arm).
#
# Generates de novo binder backbones against a fixed target, docking to a set
# of hotspot residues on the target's extracellular loops. TM bundle is left in
# the target so it acts as a steric wall (binders can't occupy membrane space).
#
# Called by run_pipeline.py, which fills the @@PLACEHOLDERS@@ from the campaign
# config. Can also be run standalone by exporting the same env vars.
#
# Required env vars (exported by the orchestrator or set manually):
#   RFD_DIR       path to RFdiffusion install
#   MODELS_DIR    path to RFdiffusion model weights
#   TARGET_PDB    fixed target PDB (full-length receptor, chain A)
#   OUT_PREFIX    output path prefix for backbones
#   CONTIG        e.g. "[A1-165/0 50-65]"  (target fixed + binder length range)
#   HOTSPOTS      e.g. "A36,A45,A47,A50,A117,A122,A128,A131"
#   NUM_DESIGNS   number of backbones to generate
#   DIFFUSER_T    diffusion timesteps (25 fast/smoketest, 50 production)
#   NOISE_SCALE   denoiser noise (0 recommended for binders -> sharper backbones)
set -eo pipefail

: "${RFD_DIR:?set RFD_DIR}"; : "${MODELS_DIR:?set MODELS_DIR}"
: "${TARGET_PDB:?set TARGET_PDB}"; : "${OUT_PREFIX:?set OUT_PREFIX}"
: "${CONTIG:?set CONTIG}"; : "${HOTSPOTS:?set HOTSPOTS}"
: "${NUM_DESIGNS:=50}"; : "${DIFFUSER_T:=50}"; : "${NOISE_SCALE:=0}"

mkdir -p "$(dirname "$OUT_PREFIX")"
cd "$RFD_DIR"

echo "[rfdiffusion] target=$TARGET_PDB contig=$CONTIG hotspots=$HOTSPOTS N=$NUM_DESIGNS T=$DIFFUSER_T"
START=$SECONDS
python scripts/run_inference.py \
  inference.output_prefix="$OUT_PREFIX" \
  inference.model_directory_path="$MODELS_DIR" \
  inference.input_pdb="$TARGET_PDB" \
  inference.num_designs="$NUM_DESIGNS" \
  "contigmap.contigs=$CONTIG" \
  "ppi.hotspot_res=[$HOTSPOTS]" \
  denoiser.noise_scale_ca="$NOISE_SCALE" \
  denoiser.noise_scale_frame="$NOISE_SCALE" \
  diffuser.T="$DIFFUSER_T"
echo "[rfdiffusion] runtime $((SECONDS-START))s"
N=$(ls "${OUT_PREFIX}"_*.pdb 2>/dev/null | wc -l)
echo "[rfdiffusion] produced $N backbones"
[ "$N" -gt 0 ] && echo "STAGE_RFD_OK"
