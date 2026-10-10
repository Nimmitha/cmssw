#!/bin/bash
# crab status for the tnp_v9 tasks that are not finished yet.
#
# First run: builds pending_tasks.txt from the current crab_* directory of every task folder
# (<stream>/<year>/<task>/crab_*; old tasks moved to <task>/removed/ are ignored).
# Every run: queries only the tasks in pending_tasks.txt. A task whose jobs are all finished is moved to
# done_tasks.txt; the others stay pending and are summarized (scheduler status, job states, exit codes).
# Full crab output of the run: status_out.txt.
#
# Usage (in tnp_v9/): ./status_pending.sh          # pending tasks only
#                     ./status_pending.sh --reset  # rebuild the list from all task folders
# After resubmitting a task into a new crab_* directory, add that directory to pending_tasks.txt (or --reset).

cd "$(dirname "$0")" || exit 1
PENDING=pending_tasks.txt
DONE=done_tasks.txt
OUT=status_out.txt

if [ "$1" == "--reset" ] || [ ! -f "$PENDING" ]; then
    ls -d muon/*/*/crab_* pdmlm/*/*/crab_* 2>/dev/null > "$PENDING"
    : > "$DONE"
    echo "pending list built: $(wc -l < "$PENDING") tasks"
fi

: > "$OUT"
still=()
n_done=0
while read -r d; do
    [ -z "$d" ] && continue
    if [ ! -d "$d" ]; then
        echo "MISSING  $d (directory not found)"
        still+=("$d")
        continue
    fi
    s=$(crab status -d "$d" 2>&1)
    { echo "Processing directory: $d"; echo "$s"; echo; } >> "$OUT"
    sched=$(echo "$s" | grep -m1 "Status on the scheduler" | awk '{print $NF}')
    jobs=$(echo "$s" | sed -n '/Jobs status:/,/^$/p' | sed 's/Jobs status://' | tr -s ' \n' ' ')
    if [ "$sched" == "COMPLETED" ] && ! echo "$jobs" | grep -qE "failed|running|idle|transferring|unsubmitted|cooloff|killed|toRetry"; then
        echo "$d" >> "$DONE"
        n_done=$((n_done + 1))
    else
        codes=$(echo "$s" | grep -oE "[0-9]+ jobs failed with exit code [0-9]+" | sed 's/ jobs failed with exit code /x/' | tr '\n' ' ')
        printf "%-55s %-10s %s %s\n" "$d" "${sched:-?}" "$jobs" "${codes:+exit codes: $codes}"
        still+=("$d")
    fi
done < "$PENDING"

printf "%s\n" "${still[@]}" | sed '/^$/d' > "$PENDING"
echo "finished this run: $n_done, still pending: $(wc -l < "$PENDING"), done in total: $(wc -l < "$DONE")"
