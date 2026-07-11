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
# Optional: NOTIFY_EMAIL=<addr> adds a single final harvest+email job.
#           CFG=<config> names the campaign config for that harvest (basename).
CFG="${CFG:-nkg7_campaign.yaml}"
REND="$WS/pipeline/rendered/$MODE"
[ -d "$REND" ] || { echo "no rendered dir $REND — run render_jobs.py first"; exit 1; }

ARMS="$*"
[ -z "$ARMS" ] && ARMS=$(ls "$REND")

JOBMAP="$WS/runs/$MODE/campaign_jobids.tsv"
mkdir -p "$(dirname "$JOBMAP")"
echo -e "arm\tstage\tjobid\tdepends_on" > "$JOBMAP"

# Collect every arm's terminal job id so a single final job can wait on them all.
TERMINAL_IDS=""
for arm in $ARMS; do
  d="$REND/$arm"
  [ -d "$d" ] || { echo "skip $arm (no rendered dir)"; continue; }
  jr=$(sbatch --parsable "$d/rfdiffusion.sbatch")
  echo -e "$arm\trfdiffusion\t$jr\t-" | tee -a "$JOBMAP"
  jm=$(sbatch --parsable --dependency=afterok:$jr "$d/proteinmpnn.sbatch")
  echo -e "$arm\tproteinmpnn\t$jm\t$jr" | tee -a "$JOBMAP"
  jb=$(sbatch --parsable --dependency=afterok:$jm "$d/boltz.sbatch")
  echo -e "$arm\tboltz\t$jb\t$jm" | tee -a "$JOBMAP"
  arm_terminal="$jb"
  # optional post-validation stages (both gated on Boltz success)
  if [ -f "$d/chai.sbatch" ]; then
    jc=$(sbatch --parsable --dependency=afterok:$jb "$d/chai.sbatch")
    echo -e "$arm\tchai\t$jc\t$jb" | tee -a "$JOBMAP"
    arm_terminal="$arm_terminal:$jc"
  fi
  if [ -f "$d/crossspecies.sbatch" ]; then
    jx=$(sbatch --parsable --dependency=afterok:$jb "$d/crossspecies.sbatch")
    echo -e "$arm\tcrossspecies\t$jx\t$jb" | tee -a "$JOBMAP"
    arm_terminal="$arm_terminal:$jx"
  fi
  TERMINAL_IDS="${TERMINAL_IDS:+$TERMINAL_IDS:}$arm_terminal"
done

# ---- Single final job: wait for ALL arms (afterany so partial failures still
# harvest), run harvest, and send ONE summary email. afterany (not afterok) so
# the notification always fires even if an arm fails. This is the only email —
# per-job mail is intentionally NOT set, to avoid dozens of messages.
if [ -n "$NOTIFY_EMAIL" ] && [ -n "$TERMINAL_IDS" ]; then
  NOTIFY="$WS/pipeline/rendered/$MODE/_notify.sbatch"
  mkdir -p "$WS/runs/$MODE/harvest"
  cat > "$NOTIFY" <<NOTIFYEOF
#!/bin/bash
#SBATCH --account=ag_iei_abdullah
#SBATCH --partition=intelsr_short
#SBATCH --time=01:00:00
#SBATCH --ntasks=1 --cpus-per-task=4 --mem=16G
#SBATCH --job-name=nkg7_${MODE}_harvest
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${NOTIFY_EMAIL}
set -eo pipefail
source /opt/software/easybuild-INTEL/software/Miniforge3/24.1.2-0/etc/profile.d/conda.sh
conda activate boltz
python $WS/pipeline/harvest.py --config $WS/pipeline/configs/$CFG \\
  --mode $MODE --ws $WS --out $WS/runs/$MODE/harvest --topn 15 \\
  > $WS/runs/$MODE/harvest/harvest.log 2>&1 || echo "harvest had errors, see log"
echo "Campaign $MODE finished. Harvest: $WS/runs/$MODE/harvest/"
cat $WS/runs/$MODE/harvest/campaign_summary.csv 2>/dev/null || true
NOTIFYEOF
  jn=$(sbatch --parsable --dependency=afterany:$TERMINAL_IDS "$NOTIFY")
  echo -e "ALL\tharvest_notify\t$jn\t$TERMINAL_IDS" | tee -a "$JOBMAP"
  echo "final harvest+email job: $jn (emails $NOTIFY_EMAIL when all arms finish)"
fi
echo "campaign submitted; job map -> $JOBMAP"
echo "watch: squeue -u \$USER"
