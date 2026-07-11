#!/bin/bash
# submit_campaign.sh — submit a full binder campaign as a chained SLURM DAG.
#
# For each arm: RFdiffusion -> (afterok) ProteinMPNN -> (afterok) Boltz-2, then
# two optional post-validation stages that both depend on Boltz:
#   Chai-1 consensus (04_chai) and cross-species re-fold (05_crossspecies).
# Uses SLURM job dependencies so each stage starts only when the previous one
# for that arm succeeds. Arms run in parallel (independent dependency chains).
# The two post-validation stages are only submitted when their sbatch scripts
# were rendered (config consensus.enabled / crossspecies.enabled).
#
# Prereq: render the sbatch scripts first:
#   python render_jobs.py --config <cfg> --mode <mode> --ws <WS> --out <WS>/pipeline/rendered
#
# Usage:
#   bash submit_campaign.sh <WS> <mode> [arm1 arm2 ...]
#   (no arm list -> all arms in the manifest)
set -eo pipefail
WS="${1:?usage: submit_campaign.sh <WS> <mode> [arms...]}"
MODE="${2:?mode: smoketest|production}"
shift 2 || true
REND="$WS/pipeline/rendered/$MODE"
[ -d "$REND" ] || { echo "no rendered dir $REND — run render_jobs.py first"; exit 1; }

ARMS="$*"
[ -z "$ARMS" ] && ARMS=$(ls "$REND")

JOBMAP="$WS/runs/$MODE/campaign_jobids.tsv"
mkdir -p "$(dirname "$JOBMAP")"
echo -e "arm\tstage\tjobid\tdepends_on" > "$JOBMAP"

for arm in $ARMS; do
  d="$REND/$arm"
  [ -d "$d" ] || { echo "skip $arm (no rendered dir)"; continue; }
  jr=$(sbatch --parsable "$d/rfdiffusion.sbatch")
  echo -e "$arm\trfdiffusion\t$jr\t-" | tee -a "$JOBMAP"
  jm=$(sbatch --parsable --dependency=afterok:$jr "$d/proteinmpnn.sbatch")
  echo -e "$arm\tproteinmpnn\t$jm\t$jr" | tee -a "$JOBMAP"
  jb=$(sbatch --parsable --dependency=afterok:$jm "$d/boltz.sbatch")
  echo -e "$arm\tboltz\t$jb\t$jm" | tee -a "$JOBMAP"
  # optional post-validation stages (both gated on Boltz success)
  if [ -f "$d/chai.sbatch" ]; then
    jc=$(sbatch --parsable --dependency=afterok:$jb "$d/chai.sbatch")
    echo -e "$arm\tchai\t$jc\t$jb" | tee -a "$JOBMAP"
  fi
  if [ -f "$d/crossspecies.sbatch" ]; then
    jx=$(sbatch --parsable --dependency=afterok:$jb "$d/crossspecies.sbatch")
    echo -e "$arm\tcrossspecies\t$jx\t$jb" | tee -a "$JOBMAP"
  fi
done
echo "campaign submitted; job map -> $JOBMAP"
echo "watch: squeue -u \$USER"
