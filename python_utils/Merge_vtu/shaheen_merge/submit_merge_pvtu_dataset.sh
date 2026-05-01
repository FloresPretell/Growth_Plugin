#!/bin/bash

set -euo pipefail

directory="/scratch/flore0a/Dataset_test2/Cases"

# Go to dataset root so dirs.txt is created/read consistently
cd "$directory"

# 1) Generate list of Simulation2D directories
python3 /project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/Postanalysis/generate_folder_list.py --path "$directory"

# 2) Count tasks
N=$(wc -l < dirs.txt)

if [[ "$N" -le 0 ]]; then
    echo "ERROR: dirs.txt is empty"
    exit 1
fi

LAST=$((N - 1))

echo "Found $N simulations -> submitting array 0-$LAST"
echo "directory=$directory"

# 3) Submit SLURM array and pass directory into the SLURM worker
sbatch --export=ALL,directory="$directory" \
    --array=0-"$LAST" \
    /project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/Merge_vtu/shaheen_merge/submit_merge_pvtu_dataset.slurm
