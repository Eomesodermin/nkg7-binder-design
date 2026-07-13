#!/bin/bash
#SBATCH --job-name=alphaflow_build2
#SBATCH --account=ag_iei_abdullah
#SBATCH --partition=mlgpu_devel
#SBATCH --gres=gpu:a40:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.out
set -o pipefail
source /opt/software/easybuild-INTEL/software/Miniforge3/24.1.2-0/etc/profile.d/conda.sh
echo "=== node: $(hostname) ==="
module load CUDA/11.8.0 && echo "CUDA: $(nvcc --version | tail -1)"
export CUDA_HOME=$EBROOTCUDA 2>/dev/null || true
conda env remove -n alphaflow -y 2>/dev/null
conda create -n alphaflow python=3.9 -y 2>&1 | tail -1
conda activate alphaflow
set -x
pip install --no-input numpy==1.21.2 pandas==1.5.3 2>&1 | tail -2
pip install --no-input torch==1.12.1+cu113 -f https://download.pytorch.org/whl/torch_stable.html 2>&1 | tail -2
pip install --no-input biopython==1.79 dm-tree==0.1.6 modelcif==0.7 ml-collections==0.1.0 scipy==1.7.1 absl-py einops 2>&1 | tail -2
# FIX: pin pytorch_lightning 1.x + torchmetrics 0.11.x (AlphaFlow-era) so the
# lightning stack does NOT drag torch up to 2.x. fair-esm/mdtraj are torch-agnostic.
pip install --no-input "pytorch_lightning==1.9.5" "torchmetrics==0.11.4" 2>&1 | tail -2
pip install --no-input fair-esm "mdtraj==1.9.9" 2>&1 | tail -2
echo "=== torch pin check BEFORE openfold ==="
python -c "import torch; print('torch', torch.__version__, 'cuda', torch.version.cuda)"
echo "=== OPENFOLD COMPILE ==="
pip install --no-input --no-build-isolation 'openfold @ git+https://github.com/aqlaboratory/openfold.git@103d037' 2>&1 | tail -20
set +x
echo "=== VERIFY ==="
python -c "import torch; print('torch', torch.__version__, 'cuda', torch.version.cuda)"
python -c "import openfold; print('openfold OK')"
python -c "import esm, mdtraj, pytorch_lightning as pl; print('deps OK, pl', pl.__version__)"
echo "ALPHAFLOW_BUILD2_DONE"
