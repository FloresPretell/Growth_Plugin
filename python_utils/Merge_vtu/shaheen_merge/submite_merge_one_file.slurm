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
#SBATCH -o /scratch/flore0a/Dataset_test2/Cases/D4_V4/Simulation2D/%x_%A_%a.out
#SBATCH -e /scratch/flore0a/Dataset_test2/Cases/D4_V4/Simulation2D/%x_%A_%a.err

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

time srun -n $_ntasks --ntasks-per-node=$_ntasks_per_node $HINT_FLAG -c $_nthreads --cpu-bind=verbose,cores   `(which pvbatch)` --force-offscreen-rendering /project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/Merge_vtu/merge_one_file.py  --root "${SIMDIR}"


# Deletion
if [[ -f "${SIMDIR}/Simulation.vtkhdf" ]]; then
    echo "SUCCESS: Simulation.vtkhdf exists."
else
    echo "ERROR: Simulation.vtkhdf was not created.VTU files will NOT be deleted."
    exit 1
fi
DELETE_LIST="${SIMDIR}/vtu_files_to_delete.txt"
if [[ ! -f "$DELETE_LIST" ]]; then
    echo "No deletion list found: $DELETE_LIST"
    echo "VTU files will NOT be deleted."
    exit 1
fi
NDELETE=$(grep -cve '^[[:space:]]*$' "$DELETE_LIST" || true)
if [[ "$NDELETE" -eq 0 ]]; then
    echo "Deletion list is empty.Nothing to delete. This usually means the VTU files were in root, not merged/."
    exit 0
fi

echo "Files marked for deletion: $NDELETE"
head "$DELETE_LIST"
echo "Deleting VTU files listed in:"
echo "$DELETE_LIST"
while IFS= read -r f; do
    [[ -z "$f" ]] && continue
    if [[ -f "$f" ]]; then
        rm -f "$f"
        echo "Deleted: $f"
    else
        echo "Already missing: $f"
    fi
done < "$DELETE_LIST"
echo "Done."