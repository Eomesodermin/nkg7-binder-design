#!/bin/bash
#SBATCH --job-name=nkg7_ga
#SBATCH --account=ag_iei_abdullah
#SBATCH --partition=mlgpu_medium
#SBATCH --gres=gpu:a40:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=20:00:00
#SBATCH --output=%x_%A_%a.out
#
# ga_refine.sbatch — partial-diffusion GENETIC ALGORITHM refinement (one arm per
# array task). Implements Tan et al. 2026 "Option I": iterate
#   select top-k backbones (global elitism) -> partial-diffuse each (mutation)
#   -> ProteinMPNN re-sequence -> Boltz fold+score -> ingest -> repeat
# for a fixed number of generations, with early-stop on no-improvement.
#
# ONE SLURM job per arm (array), throttled so only a few GPUs run at once.
# Generations are sequential within a job (gen N+1 needs gen N's selection), so
# this is the GPU-politest shape: no array of hundreds of tiny jobs.
#
# Filled by render: @@WS@@, @@ARMS@@, @@GENERATIONS@@, @@K@@, @@VARIANTS@@,
#   @@PARTIAL_T@@, @@NUM_SEQ@@, @@TOP_SEQS@@, @@DIFF_SAMPLES@@, @@PATIENCE@@
set -eo pipefail
source /opt/software/easybuild-INTEL/software/Miniforge3/24.1.2-0/etc/profile.d/conda.sh

WS=@@WS@@
ARMS=(@@ARMS@@)
ARM=${ARMS[$SLURM_ARRAY_TASK_ID]}
GENERATIONS=@@GENERATIONS@@
K=@@K@@                    # parents selected per generation (elitism breadth)
VARIANTS=@@VARIANTS@@      # partial-diffusion children per parent
PARTIAL_T=@@PARTIAL_T@@    # calibrated noise level
NUM_SEQ=@@NUM_SEQ@@        # ProteinMPNN sequences per child backbone
TOP_SEQS=@@TOP_SEQS@@      # Boltz folds only top-N seqs per backbone
DIFF_SAMPLES=@@DIFF_SAMPLES@@
PATIENCE=@@PATIENCE@@      # stop if no global-best improvement for this many gens

RFD_DIR=$WS/software/RFdiffusion
MODELS_DIR=$WS/software/RFdiffusion/models
MPNN_DIR=$WS/software/ProteinMPNN
PIPE=$WS/pipeline
GADIR=$WS/runs/refine_ga/$ARM
SEED_CSV=$WS/runs/refine_ga/refine_seed_backbones.csv
SEED_RFD=$GADIR/seed_pdbs
TARGET_SEQ_FILE=$WS/targets/$(echo $ARM | grep -q '^human' && echo human || echo mouse)_target.seq

mkdir -p $GADIR $SEED_RFD
echo "===== GA arm=$ARM gens=$GENERATIONS k=$K variants=$VARIANTS partial_T=$PARTIAL_T ====="

# --- stage this arm's seed backbone PDBs from production into one dir ---
species=$(echo $ARM | cut -d_ -f1); loop=$(echo $ARM | cut -d_ -f2)
awk -F, -v arm="$ARM" 'NR>1 && $1==arm {print $2}' $SEED_CSV | while read bb; do
  src=$WS/runs/production/$ARM/01_rfdiffusion/${bb}.pdb
  [ -f "$src" ] && cp -n "$src" $SEED_RFD/ || echo "[ga] WARN missing seed $src"
done
echo "[ga] staged $(ls $SEED_RFD/*.pdb 2>/dev/null | wc -l) seed PDBs"

# --- gen 0: init registry from production winners, select parents for gen 1 ---
conda activate boltz    # ga_select is stdlib; boltz env has python3
prev_best=$(python $PIPE/ga_select.py --mode init --arm $ARM --gen 1 --k $K \
  --registry $GADIR/registry.csv --hall $GADIR/hall.csv \
  --seed-csv $SEED_CSV --seed-rfd-dir $SEED_RFD \
  --next-parents-dir $GADIR/gen1/parents --seed-list-out $GADIR/gen1/seeds.txt \
  | sed -n 's/^GA_BEST=//p')
echo "[ga] gen0 seed best_composite=$prev_best"

stall=0
for g in $(seq 1 $GENERATIONS); do
  echo "########## GENERATION $g (arm=$ARM) ##########"
  GD=$GADIR/gen$g
  parents_dir=$GD/parents
  seeds=$(cat $GD/seeds.txt)
  [ -z "$seeds" ] && { echo "[ga] no parents -> stop"; break; }
  mkdir -p $GD/01_refine $GD/02_mpnn $GD/03_boltz

  # 1) MUTATE: partial diffusion on each selected parent backbone
  conda activate SE3nv; module load CUDA/11.7.0; export DGLBACKEND=pytorch
  RFD_DIR=$RFD_DIR MODELS_DIR=$MODELS_DIR SEED_DIR=$parents_dir \
    SEED_LIST="$seeds" OUT_PREFIX=$GD/01_refine/child \
    PARTIAL_T=$PARTIAL_T VARIANTS=$VARIANTS TARGET_LEN=165 DIFFUSER_T=50 NOISE_SCALE=0 \
    bash $PIPE/stage_refine.sh

  # 2) RE-SEQUENCE
  MPNN_DIR=$MPNN_DIR BACKBONE_DIR=$GD/01_refine OUT_DIR=$GD/02_mpnn \
    NUM_SEQ=$NUM_SEQ SAMPLING_TEMP=0.1 MODEL_NAME=v_48_020 USE_SOLUBLE=1 TARGET_CHAIN=A \
    bash $PIPE/stage_proteinmpnn.sh

  # 3) FOLD + SCORE
  conda activate boltz
  MPNN_DIR=$GD/02_mpnn TARGET_SEQ="$(cat $TARGET_SEQ_FILE)" OUT_DIR=$GD/03_boltz \
    DIFF_SAMPLES=$DIFF_SAMPLES RECYCLES=3 TOP_SEQS=$TOP_SEQS CAP=0 \
    bash $PIPE/stage_boltz.sh
  python $PIPE/filter_designs.py --boltz-dir $GD/03_boltz \
    --backbone-dir $GD/01_refine --arm $ARM --out $GD/ranked.csv

  # 4) SELECT next generation (global elitism) + improvement check
  out=$(python $PIPE/ga_select.py --mode step --arm $ARM --gen $((g+1)) --k $K \
    --registry $GADIR/registry.csv --hall $GADIR/hall.csv \
    --ranked $GD/ranked.csv --backbone-dir $GD/01_refine --prev-best $prev_best \
    --next-parents-dir $GADIR/gen$((g+1))/parents --seed-list-out $GADIR/gen$((g+1))/seeds.txt)
  echo "$out"
  new_best=$(echo "$out" | sed -n 's/^GA_BEST=//p')
  improved=$(echo "$out" | sed -n 's/^GA_IMPROVED=//p')
  echo "[ga] gen $g: best=$new_best improved=$improved (prev=$prev_best)"
  if [ "$improved" = "1" ]; then stall=0; else stall=$((stall+1)); fi
  prev_best=$new_best
  if [ "$stall" -ge "$PATIENCE" ]; then
    echo "[ga] no improvement for $PATIENCE generations -> EARLY STOP at gen $g"; break
  fi
done

echo "[ga] arm=$ARM FINAL best_composite=$prev_best"
echo "GA_ARM_DONE $ARM"
