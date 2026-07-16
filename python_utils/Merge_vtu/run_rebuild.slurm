#!/bin/bash
#SBATCH --job-name=inspect_vtkhdf
#SBATCH --partition=workq
#SBATCH --mail-user=elsanicole.florespretell@kaust.edu.sa
#SBATCH --mail-type=ALL
#SBATCH -A k10070
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --mem=80G
#SBATCH -o /scratch/flore0a/Dataset_test2/Cases/logs/%x_%j.out
#SBATCH -e /scratch/flore0a/Dataset_test2/Cases/logs/%x_%j.err

set -euo pipefail

_ntasks=1
_nthreads=1
_ntasks_per_node=1

export OMP_NUM_THREADS=$_nthreads

source /scratch/flore0a/Modules.sh

HINT_FLAG="--hint=nomultithread --mem-bind=v,none --cpu-bind=v,cores"

mkdir -p /scratch/flore0a/Dataset_test2/Cases/logs

SIMDIR="/scratch/flore0a/Dataset_test2/Cases/D8.5_V5/Simulation2D/"
SCRIPT="/project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/Merge_vtu/rebuild.py"

PVBATCH="$(which pvbatch)"

echo "SIMDIR: $SIMDIR"
echo "SCRIPT: $SCRIPT"
echo "PVBATCH: $PVBATCH"

time srun -n $_ntasks \
    --ntasks-per-node=$_ntasks_per_node \
    $HINT_FLAG \
    -c $_nthreads \
    --cpu-bind=verbose,cores \
    "$PVBATCH" \
    --force-offscreen-rendering \
    "$SCRIPT" \
    --root "$SIMDIR" \
    --make-new-temp