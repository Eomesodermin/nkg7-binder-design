#!/bin/bash
# =====================================================================
# submit_fullscale_waves.sh — wave-submission wrapper for the multi-conformer
# full-scale campaign. Submits (arm x conformer) DAG chains in bounded WAVES so
# the number of queued array tasks stays under the account's 300-submit cap,
# while keeping the throttle's worth of GPUs busy.
#
# Each (arm x conformer) chain: RFdiffusion(array) -> afterok MPNN(array) ->
# afterok Boltz(array). With coarse shards (500 bb/shard RFd=20, 500 bb/shard
# MPNN=20, 2000 folds/shard Boltz=25) each chain is 65 array tasks. The current
# config has 9 chains total (2+2+3+2 conformers across 4 arms); running all 9
# at once = 585 tasks ~2x over the 300 cap. We keep MAX_CHAINS_IN_FLIGHT chains
# (default 4 = 260 tasks, under cap, ~64 GPUs at %16) submitted, releasing the
# next as a running one clears. The per-arm consensus phase (refold->chai->cross)
# submits AFTER chains drain, so it never co-occupies the cap (peak ~106 tasks).
#
# Usage: MAX_CHAINS_IN_FLIGHT=4 bash submit_fullscale_waves.sh <WS> <REND_DIR>
# =====================================================================
set -eo pipefail
WS="${1:?need workspace}"
REND="${2:?need rendered-dir}"
MAXC="${MAX_CHAINS_IN_FLIGHT:-4}"   # 4 chains x 65 coarse-shard tasks = 260 < 300 cap, ~64 GPUs at %16
ACCT=ag_iei_abdullah
JOBMAP="$WS/runs/fullscale/chain_jobids.tsv"
mkdir -p "$(dirname "$JOBMAP")"; : > "$JOBMAP"

# Discover sharded chains = subdirs named <arm>__<confid> (EXCLUDE <arm>__consensus,
# which holds chai/crossspecies and is submitted separately, held on its chains).
CHAINS=($(ls -d "$REND"/*/ 2>/dev/null | xargs -n1 basename | grep -v '__consensus$'))
echo "[waves] ${#CHAINS[@]} sharded chains to run, max $MAXC in flight, cap-aware"

# count our currently-submitted array tasks (PENDING+RUNNING), to stay < 300
submitted_tasks() {
  squeue -u "$USER" -h -t PENDING,RUNNING -r 2>/dev/null | wc -l
}

launch_chain() {
  local d="$REND/$1"
  local jr jm jb
  jr=$(sbatch --parsable "$d/rfdiffusion.sbatch")
  jm=$(sbatch --parsable --dependency=afterok:$jr "$d/proteinmpnn.sbatch")
  jb=$(sbatch --parsable --dependency=afterok:$jm "$d/boltz.sbatch")
  echo -e "$1\t$jr\t$jm\t$jb" >> "$JOBMAP"
  echo "[waves] launched $1 : rfd=$jr mpnn=$jm boltz=$jb"
}

idx=0; inflight=()
while [ $idx -lt ${#CHAINS[@]} ] || [ ${#inflight[@]} -gt 0 ]; do
  # prune finished chains from inflight (boltz job no longer in queue)
  new_inflight=()
  for entry in "${inflight[@]}"; do
    cbb="${entry##*:}"
    if squeue -j "$cbb" -h 2>/dev/null | grep -q .; then new_inflight+=("$entry"); fi
  done
  inflight=("${new_inflight[@]}")

  # launch while we have room in chain-slots AND task-cap headroom
  while [ ${#inflight[@]} -lt $MAXC ] && [ $idx -lt ${#CHAINS[@]} ]; do
    # cap guard: don't submit a ~102-task chain if it would breach ~290
    if [ "$(submitted_tasks)" -gt 260 ]; then break; fi
    launch_chain "${CHAINS[$idx]}"
    boltz_id=$(tail -1 "$JOBMAP" | cut -f4)
    inflight+=("${CHAINS[$idx]}:$boltz_id")
    idx=$((idx+1))
  done
  [ $idx -lt ${#CHAINS[@]} ] || [ ${#inflight[@]} -gt 0 ] && sleep 120
done

echo "[waves] all ${#CHAINS[@]} sharded chains submitted; job map: $JOBMAP"

# --- Consensus phase: per arm, a 5-stage chain held on afterok of ALL that arm's
# conformer-chain boltz jobs (pools designs across conformers, runs once):
#   refold_prep -> refold_fold(array) -> refold_rank -> chai -> crossspecies
echo "[waves] submitting per-arm consensus (ds=5 refold -> chai -> cross), held on chains..."
LAST_CROSS_IDS=()
for cdir in $(ls -d "$REND"/*__consensus/ 2>/dev/null | xargs -n1 basename); do
  arm="${cdir%%__consensus}"
  deps=$(awk -F'\t' -v a="${arm}__" '$1 ~ "^"a {print $4}' "$JOBMAP" | paste -sd: -)
  if [ -z "$deps" ]; then echo "[waves] WARN no chains for $arm, skipping consensus"; continue; fi
  jp=$(sbatch --parsable --dependency=afterok:$deps "$REND/$cdir/refold_prep.sbatch")
  jf=$(sbatch --parsable --dependency=afterok:$jp "$REND/$cdir/refold_fold.sbatch")
  jr=$(sbatch --parsable --dependency=afterok:$jf "$REND/$cdir/refold_rank.sbatch")
  jc=$(sbatch --parsable --dependency=afterok:$jr "$REND/$cdir/chai.sbatch")
  jx=$(sbatch --parsable --dependency=afterok:$jc "$REND/$cdir/crossspecies.sbatch")
  echo -e "${arm}__consensus\trefprep=$jp\trefold=$jf\trefrank=$jr\tchai=$jc\tcross=$jx" >> "$JOBMAP"
  echo "[waves] $arm consensus: refprep=$jp refold=$jf refrank=$jr chai=$jc cross=$jx (afterok:$deps)"
  LAST_CROSS_IDS+=("$jx")
done

# --- Email-on-completion: one tiny job held afterany ALL cross jobs, emails the user.
if [ ${#LAST_CROSS_IDS[@]} -gt 0 ] && [ -n "${NOTIFY_EMAIL:-}" ]; then
  alldeps=$(IFS=:; echo "${LAST_CROSS_IDS[*]}")
  jn=$(sbatch --parsable --dependency=afterany:$alldeps \
        --account=ag_iei_abdullah --partition=intelsr_short --time=00:05:00 \
        --ntasks=1 --cpus-per-task=1 --mem=2G --job-name=nkg7_done \
        --mail-user="$NOTIFY_EMAIL" --mail-type=END,FAIL \
        --wrap="echo 'NKG7 full-scale campaign complete. Results under $WS/runs/fullscale/*/{04_chai,05_crossspecies,06_refold}.'")
  echo -e "notify\t-\t-\t-\t-\temail=$jn" >> "$JOBMAP"
  echo "[waves] email-on-completion job=$jn (afterany:$alldeps) -> $NOTIFY_EMAIL"
elif [ -z "${NOTIFY_EMAIL:-}" ]; then
  echo "[waves] NOTE: NOTIFY_EMAIL unset -> no completion email submitted"
fi
echo "[waves] complete; job map: $JOBMAP"
