# Conda environments (marvin)

Four conda environments back the pipeline. They mix conda + pip heavily (PyTorch
CUDA wheels, Boltz-2, Chai-1 installed via pip), so the `*.pip-freeze.txt` files
are the load-bearing record; the `*.environment.yml` files capture the conda
layer. Rebuild an env with the pip-freeze after creating a base conda env at the
right Python version.

| Env | Python | Purpose | Key packages |
|---|---|---|---|
| `SE3nv` | 3.9 | RFdiffusion (design + partial-diffusion refinement) | rfdiffusion 1.1.0, torch 1.12.1+cu116, dgl 1.0.2+cu116, numpy 1.26.4 |
| `boltz` | 3.11 | Boltz-2 validation (folding/scoring) | boltz 2.2.1, torch 2.13 |
| `chai` | 3.11 | Chai-1 orthogonal consensus | chai_lab 0.6.1, torch 2.6 |
| `openmm-env` | 3.11 | loop-flexibility MD ensemble | openmm 8.2, pdbfixer |

**marvin AVX-512 gotcha:** the easybuild-INTEL module stack is AVX-512-compiled and
SIGILLs (exit 132) on the AMD mlgpu (A40) nodes, which have AVX2 but not AVX-512.
That is why Boltz-2/Chai are pip-installed (AVX2-safe wheels) rather than loaded
from the ColabFold/AlphaFold/Boltz-1 modules.

**SE3nv runtime notes:** needs `module load CUDA/11.7.0` (for libcusparse.so.11) and
`export DGLBACKEND=pytorch` before running.
