#!/bin/bash

SIMDIR="/Users/flore0a/UG4-IA/ug4/apps/Neuronal_Process/Scripts_slurm_sh/Simulation2D"
PYTHON=/opt/anaconda3/envs/paraview_env/bin/python

"$PYTHON" merge_one_file.py  --root "${SIMDIR}"


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