# Conda environments

Four conda environments back the pipeline. They mix conda and pip heavily (PyTorch CUDA wheels,
Boltz-2 and Chai-1 installed with pip), so the `*.pip-freeze.txt` files are the load-bearing
record and the `*.environment.yml` files capture the conda layer. To rebuild one, create a base
environment at the right Python version and install from the pip-freeze file.

| Env | Python | Purpose | Key packages |
|---|---|---|---|
| `SE3nv` | 3.9 | RFdiffusion (design and partial-diffusion refinement) and ProteinMPNN | rfdiffusion 1.1.0, torch 1.12.1+cu116, dgl 1.0.2+cu116, numpy 1.26.4 |
| `boltz` | 3.11 | Boltz-2 validation (folding and scoring) | boltz 2.2.1, torch 2.13 |
| `chai` | 3.11 | Chai-1 orthogonal consensus | chai_lab 0.6.1, torch 2.6 |
| `openmm-env` | 3.11 | Loop-flexibility MD ensembles | openmm 8.2, pdbfixer |

Entries that pointed at conda-build directories on the original machine (`@ file://...`) were
removed from the pip-freeze files. Those packages are covered by the conda layer.

**AVX-512 note.** The cluster's EasyBuild module stack is compiled for AVX-512 and crashes
(SIGILL, exit 132) on its AMD GPU nodes, which have AVX2 only. That is why Boltz-2 and Chai-1 are
pip-installed with generic wheels instead of loaded from the ColabFold, AlphaFold or Boltz-1
modules. If your nodes support AVX-512 the modules may work.

**SE3nv runtime notes.** It needs `module load CUDA/11.7.0` (for `libcusparse.so.11`) and
`export DGLBACKEND=pytorch` before running.
