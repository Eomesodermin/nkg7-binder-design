#!/bin/bash
#SBATCH --job-name=nkg7_refine_calib
#SBATCH --account=ag_iei_abdullah
#SBATCH --partition=mlgpu_short
#SBATCH --gres=gpu:a40:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%j.out
set -eo pipefail
source /opt/software/easybuild-INTEL/software/Miniforge3/24.1.2-0/etc/profile.d/conda.sh

WS=/lustre/scratch/data/dcorvino_hpc-nkg7_binder_design
RFD_DIR=$WS/software/RFdiffusion
MODELS_DIR=$WS/software/RFdiffusion/models
MPNN_DIR=$WS/software/ProteinMPNN
PIPE=$WS/pipeline
CALIB=$WS/runs/refine_calib
mkdir -p $CALIB/seeds

# stage the two calibration seed PDBs into one dir
cp $WS/runs/production/mouse_ecl1/01_rfdiffusion/mouse_ecl1_355.pdb $CALIB/seeds/
cp $WS/runs/production/human_ecl2/01_rfdiffusion/human_ecl2_482.pdb $CALIB/seeds/

# (seed, species target seq) pairs
declare -A TGT
TGT[mouse_ecl1_355]=$WS/targets/mouse_target.seq
TGT[human_ecl2_482]=$WS/targets/human_target.seq

for PT in 5 10 20; do
  for seed in mouse_ecl1_355 human_ecl2_482; do
    d=$CALIB/${seed}_pt${PT}
    mkdir -p $d/01_refine $d/02_mpnn $d/03_boltz
    echo "=================== $seed partial_T=$PT ==================="

    # --- refine (RFdiffusion partial diffusion) ---
    conda activate SE3nv; module load CUDA/11.7.0; export DGLBACKEND=pytorch
    RFD_DIR=$RFD_DIR MODELS_DIR=$MODELS_DIR SEED_DIR=$CALIB/seeds \
      SEED_LIST="$seed" OUT_PREFIX=$d/01_refine/${seed} \
      PARTIAL_T=$PT VARIANTS=10 TARGET_LEN=165 DIFFUSER_T=50 NOISE_SCALE=0 \
      bash $PIPE/stage_refine.sh

    # --- ProteinMPNN re-sequence ---
    MPNN_DIR=$MPNN_DIR BACKBONE_DIR=$d/01_refine OUT_DIR=$d/02_mpnn \
      NUM_SEQ=8 SAMPLING_TEMP=0.1 MODEL_NAME=v_48_020 USE_SOLUBLE=1 TARGET_CHAIN=A \
      bash $PIPE/stage_proteinmpnn.sh

    # --- Boltz fold (top-2 seqs/backbone, ds=1) ---
    conda activate boltz
    MPNN_DIR=$d/02_mpnn TARGET_SEQ="$(cat ${TGT[$seed]})" OUT_DIR=$d/03_boltz \
      DIFF_SAMPLES=1 RECYCLES=3 TOP_SEQS=2 CAP=0 \
      bash $PIPE/stage_boltz.sh

    # --- score ---
    python $PIPE/filter_designs.py --boltz-dir $d/03_boltz \
      --backbone-dir $d/01_refine --arm ${seed}_pt${PT} \
      --out $CALIB/${seed}_pt${PT}_ranked.csv || echo "filter failed $seed pt$PT"
  done
done
echo "REFINE_CALIB_DONE"
