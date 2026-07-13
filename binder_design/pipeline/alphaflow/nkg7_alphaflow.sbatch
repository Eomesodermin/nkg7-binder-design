#!/bin/bash
#SBATCH --job-name=nkg7_alphaflow
#SBATCH --account=ag_iei_abdullah
#SBATCH --partition=mlgpu_short
#SBATCH --gres=gpu:a40:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%j.out
#
# nkg7_alphaflow.sbatch — Option-2 DL-ensemble cross-check. MSA-backed AlphaFlow
# MD+Templates (12l-distilled), 25 samples per target = an ensemble orthogonal to
# the explicit-solvent MD. Both NKG7 species in one job (2 short seqs).
set -eo pipefail
source /opt/software/easybuild-INTEL/software/Miniforge3/24.1.2-0/etc/profile.d/conda.sh
module load CUDA/11.8.0
conda activate alphaflow
AF=/lustre/scratch/data/dcorvino_hpc-nkg7_binder_design/software/alphaflow
RUN=/lustre/scratch/data/dcorvino_hpc-nkg7_binder_design/runs/alphaflow
cd $AF
echo "===== AlphaFlow MSA-backed on $(hostname) ====="
python predict.py --mode alphafold   --input_csv $RUN/input/nkg7.csv   --msa_dir $RUN/msa   --templates_dir $RUN/templates   --weights $AF/weights/alphaflow_12l_md_templates_distilled_202406.pt   --samples 25 --steps 10   --outpdb $RUN/outpdb
echo "=== output conformers ==="
ls $RUN/outpdb/ | head
echo "ALPHAFLOW_INFER_DONE"
