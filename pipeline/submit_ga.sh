#!/bin/bash
# submit_ga.sh — render + submit the partial-diffusion GA refinement campaign.
#
# Usage: NOTIFY_EMAIL=you@uni-bonn.de bash submit_ga.sh <WS>
#
# Submits ONE array job (one task per arm), throttled %2 (<=2 GPUs at once), plus
# a single dependent harvest+notify job that emails ONCE when all arms finish.
# Good-cluster-citizen shape: 1 array job + 1 notify job = 2 submissions total,
# no per-generation job spam, throttled GPU use, early-stop frees GPUs.
set -eo pipefail
WS=${1:?usage: submit_ga.sh <WS>}
: "${NOTIFY_EMAIL:?set NOTIFY_EMAIL for completion notification}"

# --- GA hyperparameters (calibrated) ---
ARMS_LIST="human_ecl1 human_ecl2 mouse_ecl1 mouse_ecl2"   # composite arms dropped
GENERATIONS=12
K=4                 # parents/generation (global elitism breadth)
VARIANTS=8          # partial-diffusion children per parent
PARTIAL_T=10        # calibrated: best fold quality without pt20 drift
NUM_SEQ=8
TOP_SEQS=2          # Boltz folds top-2 seqs/child backbone
DIFF_SAMPLES=1
PATIENCE=3          # early-stop after 3 gens with no global-best improvement
THROTTLE=2          # max concurrent array tasks (<=2 GPUs at once)

PIPE=$WS/pipeline
RUNDIR=$WS/runs/refine_ga
mkdir -p $RUNDIR

# arms as a bash-array literal for the sbatch script
narms=$(echo $ARMS_LIST | wc -w)
armquoted=$(for a in $ARMS_LIST; do echo -n "\"$a\" "; done)

# --- render the array driver from the template ---
sed -e "s|@@WS@@|$WS|g" \
    -e "s|@@ARMS@@|$armquoted|g" \
    -e "s|@@GENERATIONS@@|$GENERATIONS|g" \
    -e "s|@@K@@|$K|g" \
    -e "s|@@VARIANTS@@|$VARIANTS|g" \
    -e "s|@@PARTIAL_T@@|$PARTIAL_T|g" \
    -e "s|@@NUM_SEQ@@|$NUM_SEQ|g" \
    -e "s|@@TOP_SEQS@@|$TOP_SEQS|g" \
    -e "s|@@DIFF_SAMPLES@@|$DIFF_SAMPLES|g" \
    -e "s|@@PATIENCE@@|$PATIENCE|g" \
    $PIPE/ga_refine.sbatch > $RUNDIR/ga_refine.rendered.sbatch
bash -n $RUNDIR/ga_refine.rendered.sbatch

# --- submit the arm array, throttled ---
last=$((narms-1))
GA_JID=$(sbatch --parsable --array=0-${last}%${THROTTLE} \
  --chdir=$RUNDIR $RUNDIR/ga_refine.rendered.sbatch)
echo "[submit] GA array job: $GA_JID  (arms: $ARMS_LIST, throttle %$THROTTLE)"

# --- harvest + single email, depends on ALL arms finishing (afterany) ---
NOTIFY_JID=$(sbatch --parsable --dependency=afterany:${GA_JID} \
  --account=ag_iei_abdullah --partition=intelsr_short --time=00:20:00 \
  --job-name=nkg7_ga_notify --mail-type=END,FAIL --mail-user=$NOTIFY_EMAIL \
  --output=$RUNDIR/ga_notify_%j.out \
  --wrap="source /opt/software/easybuild-INTEL/software/Miniforge3/24.1.2-0/etc/profile.d/conda.sh; conda activate boltz; python $PIPE/ga_harvest.py --ga-dir $RUNDIR --arms '$ARMS_LIST' --out $RUNDIR/ga_summary.csv")
echo "[submit] harvest+notify job: $NOTIFY_JID (afterany:$GA_JID, emails $NOTIFY_EMAIL)"

printf "arm_array\t%s\nnotify\t%s\n" "$GA_JID" "$NOTIFY_JID" > $RUNDIR/ga_jobids.tsv
echo "[submit] job map -> $RUNDIR/ga_jobids.tsv"
