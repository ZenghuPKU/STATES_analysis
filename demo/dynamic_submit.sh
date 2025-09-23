#!/bin/bash

task_dir="local_position_subtile_tasks"
scripts=(${task_dir}/local_task_Position*_subtile*.sh)

if [[ ${#scripts[@]} -eq 0 ]]; then
    echo "Error: No task scripts found in ${task_dir}. Exiting."
    exit 1
fi

max_running_jobs=25
total_scripts=${#scripts[@]}

current_index=0

submit_task() {
    local task_script=$1
    sbatch "$task_script" > /dev/null 2>&1
    if [[ $? -eq 0 ]]; then
        echo "Submitted: $task_script"
    else
        echo "Error: Failed to submit $task_script"
    fi
}

while [[ $current_index -lt $total_scripts ]]; do
    running_tasks=$(squeue -u $USER | grep -E "fat4way" | wc -l)

    if [[ $running_tasks -lt $max_running_jobs ]]; then
        submit_task "${scripts[current_index]}"
        current_index=$((current_index + 1))
    else
        echo "Currently running $running_tasks tasks. Waiting for a slot..."
        sleep 10
    fi
done

echo "All tasks submitted. Waiting for completion..."
while true; do
    running_tasks=$(squeue -u $USER | grep -E "fat4way" | wc -l)
    if [[ $running_tasks -eq 0 ]]; then
        echo "All tasks completed!"
        break
    else
        echo "$running_tasks tasks are still running. Checking again in 10 seconds..."
        sleep 10
    fi
done
