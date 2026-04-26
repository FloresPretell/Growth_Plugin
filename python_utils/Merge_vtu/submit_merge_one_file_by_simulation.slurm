#!/bin/bash
#SBATCH --job-name=MergeToVTKHDF_by_Simulation
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
SIMDIR="${ROOT%/}" # here is all the outputs (Simulations2D)

time srun -n $_ntasks --ntasks-per-node=$_ntasks_per_node $HINT_FLAG -c $_nthreads --cpu-bind=verbose,cores   `(which pvbatch)` --force-offscreen-rendering /project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/Merge_vtu/merge_one_file.py  --root "${SIMDIR}"

if [ -f "${SIMDIR}/simulation.vtkhdf" ] || [ -f "${SIMDIR}/simulation.xmf" ]; then
  rm -f "${SIMDIR}/merged"/Output_t*.vtu
fi