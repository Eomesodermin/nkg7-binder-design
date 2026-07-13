#!/bin/bash
#SBATCH --job-name=nkg7_mdexpl
#SBATCH --account=ag_iei_abdullah
#SBATCH --partition=mlgpu_medium
#SBATCH --gres=gpu:a40:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=20:00:00
#SBATCH --output=%x_%j.out
#SBATCH --array=0-1
#
# md_explicit.sbatch — explicit-solvent (TIP3P/PME/NPT) loop-ensemble MD, the
# Option-3 escalation. One array task per species (0=human, 1=mouse), each a
# single-GPU 50 ns run. Independent of the DL-ensemble (AlphaFlow) job.
set -eo pipefail
source /opt/software/easybuild-INTEL/software/Miniforge3/24.1.2-0/etc/profile.d/conda.sh
conda activate openmm-env

WS=/lustre/scratch/data/dcorvino_hpc-nkg7_binder_design
SPECIES=(human mouse)
sp=${SPECIES[$SLURM_ARRAY_TASK_ID]}
declare -A PDB
PDB[human]=$WS/targets/NKG7_human_Q16617_AF.pdb
PDB[mouse]=$WS/targets/NKG7_mouse_Q99PA5_AF.pdb

echo "===== explicit-solvent MD: $sp on $(hostname) ====="
python $WS/pipeline/make_loop_ensemble_explicit.py \
  --pdb ${PDB[$sp]} --loops 30-60 113-132 \
  --n-conformers 25 --production-ns 50 --equil-ns 0.5 \
  --seed 2024 --out $WS/runs/md_explicit/$sp
echo "MD_EXPLICIT_DONE $sp"
