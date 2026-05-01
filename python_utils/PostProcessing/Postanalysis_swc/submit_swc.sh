#!/bin/bash

set -euo pipefail
PYTHON_BIN="/scratch/flore0a/iops/miniconda3/envs/swc_env/bin/python"
directorio=/scratch/flore0a/Dataset_test2/Cases

# Stage 1: generate folder list
"${PYTHON_BIN}" generate_folder_list_swc.py --path "${directorio}"

# Count tasks
N=$(wc -l < dirs.txt)
if [[ "$N" -le 0 ]]; then
  echo "ERROR: dirs.txt está vacío"
  exit 1
fi
LAST=$((N-1))

echo "Found $N simulations -> submitting SWC export array 0-$LAST"

# Stage 1: submit SWC export array, capture job ID
SWC_JOB=$(sbatch --array=0-"$LAST" Submit_export_swc.slurm | awk '{print $NF}')
echo "SWC export job: ${SWC_JOB}"

# Stage 2: submit morphometrics job, runs only after all array tasks complete
MET_JOB=$(sbatch --dependency=afterok:"${SWC_JOB}" Submit_morphometrics.slurm | awk '{print $NF}')
echo "Morphometrics job: ${MET_JOB}  (depends on ${SWC_JOB})"
