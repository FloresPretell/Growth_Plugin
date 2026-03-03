#!/bin/bash
#SBATCH --job-name=Export_png
#SBATCH --partition=workq
#SBATCH --mail-user=elsanicole.florespretell@kaust.edu.sa
#SBATCH --mail-type=ALL
#SBATCH -A k10070
#SBATCH --nodes=2
#SBATCH --time=23:59:00
#SBATCH --array=0-5
#SBATCH --mem=200G
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail
mkdir -p logs

# set script vars
_ntasks=8
_nthreads=48
_ntasks_per_node=5

#OpenMP settings:
export OMP_NUM_THREADS=$_nthreads

# load the paraview environment
module load paraview
source /scratch/flore0a/Modules.sh

HINT_FLAG="--hint=nomultithread  --mem-bind=v,none --cpu-bind=v,cores"

# Read dirs.txt
SIMDIR=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" dirs.txt)

echo "Task ${SLURM_ARRAY_TASK_ID} exporting ${SIMDIR}"

time srun -n $_ntasks --ntasks-per-node=$_ntasks_per_node $HINT_FLAG -c $_nthreads --cpu-bind=verbose,cores   `(which pvbatch)` --force-offscreen-rendering ./main_export.py  --simdir "${SIMDIR}" --variable_input "t"
