#!/bin/bash
#SBATCH --job-name=merge_pvtu_to_vtu_individual
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
SIMDIR="/scratch/flore0a/Dataset_test4/Cases/D10_V0p0075/Simulation2D/"

time srun -n $_ntasks --ntasks-per-node=$_ntasks_per_node $HINT_FLAG -c $_nthreads --cpu-bind=verbose,cores   `(which pvbatch)` --force-offscreen-rendering /project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/Merge_vtu/merge_pvtu_to_vtu.py  --root "${SIMDIR}" --list_to_delete "${SIMDIR}/files_to_delete.txt"

### STEP 2
# delete fragments listed by python
SIMDIR=${SIMDIR%/} # remove trailing slash if any
LIST="${SIMDIR}/files_to_delete.txt"

if [[ -s "$LIST" ]]; then
  while IFS= read -r pvtu; do
    [[ -z "$pvtu" ]] && continue

    t_str="$(basename "$pvtu")"
    t_str="${t_str#Output_t}"
    t_str="${t_str%.pvtu}"

    merged_vtu="${SIMDIR}/merged/Output_t${t_str}.vtu"
    if [[ ! -s "$merged_vtu" ]]; then
      echo "[delete] SKIP missing or empty merged file: $merged_vtu"
      continue
    fi

    echo "[delete] deleting ${pvtu}"
    rm -f -- "$pvtu"

    echo "[delete] deleting partition files Output_p*_t${t_str}.vtu"
    find "$SIMDIR" -maxdepth 1 -type f -name "Output_p*_t${t_str}.vtu" -print -delete

  done < "$LIST"

  rm -f -- "$LIST"
else
  echo "[delete] WARNING: list file missing or empty"
fi