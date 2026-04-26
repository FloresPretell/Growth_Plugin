#!/bin/bash
#SBATCH --job-name=merge_to_vtkhdf
#SBATCH --partition=workq
#SBATCH --mail-user=elsanicole.florespretell@kaust.edu.sa
#SBATCH --mail-type=ALL
#SBATCH -A k10070
#SBATCH --nodes=1
#SBATCH --time=23:59:00
#SBATCH --array=0-0
#SBATCH --mem=200G
#SBATCH -o /dev/null
#SBATCH -e /dev/null

set -euo pipefail

# set script vars
_ntasks=1
_nthreads=1
_ntasks_per_node=1

#OpenMP settings:
export OMP_NUM_THREADS=$_nthreads

# load the paraview environment
source /scratch/flore0a/Modules.sh
HINT_FLAG="--hint=nomultithread  --mem-bind=v,none --cpu-bind=v,cores"
SIMDIR="/scratch/flore0a/Dataset_test2/Cases/D4_V4/Simulation2D/"

time srun -n $_ntasks --ntasks-per-node=$_ntasks_per_node $HINT_FLAG -c $_nthreads --cpu-bind=verbose,cores   `(which pvbatch)` --force-offscreen-rendering merge_one_file.py  --root "${SIMDIR}"