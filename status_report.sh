#!/bin/bash
# Where does every case stand, and what must be (re)submitted? Works on the login node.
#   bash status_report.sh
# Categories:
#   RUNNING   in the queue now
#   DONE      finished with the current code (has the lower-trench metric c_deep)
#   DONE-OLD  finished before the comparison update (no c_deep): rerun with --force for the comparison
#   STOPPED   checkpoint left behind (time limit, maintenance, crash): resume it
#   FAILED    folder exists but no result and no checkpoint: submit again
#   NOT-RUN   never started
cd "$(dirname "$0")"

declare -A group
order=(RUNNING DONE DONE-OLD STOPPED FAILED NOT-RUN)
queued=$(squeue -u "$USER" -h -o %j 2>/dev/null | sort -u)

while IFS=, read -r case _ _ _ _ _ _ _ _ walltime _; do
    dir="results/$case"
    if grep -qx "$case" <<<"$queued"; then
        st=RUNNING
    elif [[ -f "$dir/summary.json" ]]; then
        # A0 (original mesh) and T (wafer-scale tank) have no deep-trench metric by design
        if head -1 "$dir/history.csv" 2>/dev/null | grep -q "c_deep" || [[ "$case" == A0_* || "$case" == T_* ]]; then
            st=DONE
        else
            st=DONE-OLD
        fi
    elif [[ -f "$dir/checkpoint.json" ]]; then
        st=STOPPED
    elif [[ -d "$dir" ]]; then
        st=FAILED
    else
        st=NOT-RUN
    fi
    group[$st]+="$case@$walltime "
done < <(tail -n +2 cases.csv)

for st in "${order[@]}"; do
    n=$(wc -w <<<"${group[$st]}")
    printf "%-9s %3d  %s\n" "$st" "$n" "$(tr ' ' '\n' <<<"${group[$st]}" | cut -d@ -f1 | head -6 | paste -sd' ')$( ((n > 6)) && echo ' ...')"
done

echo
echo "---- commands (run from ~/3d_nand_2) ----"
regex() { tr ' ' '\n' <<<"$1" | cut -d@ -f1 | grep . | paste -sd'|'; }
if [[ -n "${group[STOPPED]}" ]]; then
    echo "# resume stopped runs from their checkpoints:"
    for cw in ${group[STOPPED]}; do
        echo "sbatch --job-name=${cw%%@*} --time=${cw##*@} run_case.slurm ${cw%%@*} --resume"
    done
fi
if [[ -n "${group[FAILED]}${group[NOT-RUN]}" ]]; then
    echo "# start failed / never-run cases:"
    echo "bash submit_cases.sh '^($(regex "${group[FAILED]} ${group[NOT-RUN]}"))\$'"
fi
if [[ -n "${group[DONE-OLD]}" ]]; then
    echo "# rerun cases finished with the old code (needed for the deep-trench comparison):"
    echo "bash submit_cases.sh '^($(regex "${group[DONE-OLD]}"))\$' --force"
fi
[[ -z "${group[STOPPED]}${group[FAILED]}${group[NOT-RUN]}${group[DONE-OLD]}${group[RUNNING]}" ]] && \
    echo "# everything is DONE: run analyze.py on a compute node"
