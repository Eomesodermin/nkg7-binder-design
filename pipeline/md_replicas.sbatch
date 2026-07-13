#!/bin/bash
#SBATCH --job-name=nkg7_mdrep
#SBATCH --account=ag_iei_abdullah
#SBATCH --partition=mlgpu_short
#SBATCH --gres=gpu:a40:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%j.out
#
# md_replicas.sbatch — basin-convergence robustness check for the NKG7 loop
# ensemble. Runs N independent restrained-MD trajectories (different random
# seeds) per species, each longer than the original diagnostic run, then
# harvests conformers from every replica. A downstream analysis compares
# per-replica loop RMSF and cross-replica structural spread to decide whether
# the loops sample ONE basin (our 12-conformer set is representative -> proceed)
# or MULTIPLE basins (escalate to DL-ensemble + explicit-solvent MD).
set -eo pipefail
source /opt/software/easybuild-INTEL/software/Miniforge3/24.1.2-0/etc/profile.d/conda.sh
conda activate openmm-env

WS=/lustre/scratch/data/dcorvino_hpc-nkg7_binder_design
PIPE=$WS/pipeline
OUT=$WS/runs/md_replicas
mkdir -p $OUT

# Longer sampling than the diagnostic run: 100 ps equil, snapshots every 50 ps
# (less correlated), 20 conformers/replica = ~1 ns production per replica.
NREP=4
NCONF=20
EQUIL=50000        # 100 ps
INTERVAL=25000     # 50 ps between snapshots

declare -A LOOPS
LOOPS[human]="30-60 113-132"
LOOPS[mouse]="30-60 113-132"
declare -A PDB
PDB[human]=$WS/targets/NKG7_human_Q16617_AF.pdb
PDB[mouse]=$WS/targets/NKG7_mouse_Q99PA5_AF.pdb

for sp in human mouse; do
  for rep in $(seq 1 $NREP); do
    seed=$((rep*7919))   # distinct prime-spaced seeds
    d=$OUT/${sp}_rep${rep}
    echo "===== $sp replica $rep (seed=$seed) ====="
    python $PIPE/make_loop_ensemble.py \
      --pdb ${PDB[$sp]} --loops ${LOOPS[$sp]} \
      --n-conformers $NCONF --equil-steps $EQUIL --sample-interval $INTERVAL \
      --seed $seed --out $d
  done
done
echo "MD_REPLICAS_DONE"
