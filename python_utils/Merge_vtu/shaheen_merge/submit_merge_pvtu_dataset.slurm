#!/bin/bash
#SBATCH --job-name=merge_pvtu_to_vtu
#SBATCH --partition=workq
#SBATCH --mail-user=elsanicole.florespretell@kaust.edu.sa
#SBATCH --mail-type=ALL
#SBATCH -A k10070
#SBATCH --nodes=1
#SBATCH --time=23:59:00
#SBATCH --mem=200G
#SBATCH -o /scratch/flore0a/Dataset_test4/Cases/logs/%x_%A_%a.out
#SBATCH -e /scratch/flore0a/Dataset_test4/Cases/logs/%x_%A_%a.err

set -euo pipefail

# ----------------------------
# Basic task settings
# ----------------------------
_ntasks=1
_nthreads=1
_ntasks_per_node=1

export OMP_NUM_THREADS=$_nthreads

# ----------------------------
# Load ParaView environment
# ----------------------------
source /scratch/flore0a/Modules.sh

HINT_FLAG="--hint=nomultithread --mem-bind=v,none --cpu-bind=v,cores"

# ----------------------------
# Dataset root passed from submitter
# ----------------------------
CASES_ROOT="${directory:?ERROR: directory was not passed from the submit script}"
DIRS_FILE="${CASES_ROOT}/dirs.txt"

mkdir -p "${CASES_ROOT}/logs"

if [[ ! -f "$DIRS_FILE" ]]; then
    echo "ERROR: dirs.txt not found:"
    echo "$DIRS_FILE"
    exit 1
fi

# ----------------------------
# Select Simulation2D directory for this array task
# ----------------------------
CASE_DIR=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$DIRS_FILE")

if [[ -z "${CASE_DIR:-}" ]]; then
    echo "ERROR: No directory found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# Your dirs.txt already contains Simulation2D directories
SIMDIR="${CASE_DIR%/}/"

echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
echo "CASES_ROOT: $CASES_ROOT"
echo "CASE_DIR: $CASE_DIR"
echo "SIMDIR: $SIMDIR"

if [[ ! -d "$SIMDIR" ]]; then
    echo "ERROR: Simulation2D directory does not exist:"
    echo "$SIMDIR"
    exit 1
fi

cd "$SIMDIR"

PVBATCH="$(which pvbatch)"
echo "Using pvbatch:"
echo "$PVBATCH"

# ----------------------------
# STEP 1: merge pvtu fragments into vtu files
# ----------------------------
time srun -n $_ntasks \
    --ntasks-per-node=$_ntasks_per_node \
    $HINT_FLAG \
    -c $_nthreads \
    --cpu-bind=verbose,cores \
    "$PVBATCH" \
    --force-offscreen-rendering \
    /project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/Merge_vtu/merge_pvtu_to_vtu.py \
    --root "$SIMDIR" \
    --list_to_delete "${SIMDIR}/files_to_delete.txt"

# ----------------------------
# STEP 2: delete fragments listed by python
# ----------------------------
SIMDIR="${SIMDIR%/}"
LIST="${SIMDIR}/files_to_delete.txt"

if [[ -s "$LIST" ]]; then
    while IFS= read -r pvtu; do
        [[ -z "$pvtu" ]] && continue

        t_str="$(basename "$pvtu")"
        t_str="${t_str#Output_t}"
        t_str="${t_str%.pvtu}"

        merged_vtu="${SIMDIR}/merged/Output_t${t_str}.vtu"

        if [[ ! -s "$merged_vtu" ]]; then
            echo "[delete] SKIP missing or empty merged file:"
            echo "$merged_vtu"
            continue
        fi

        echo "[delete] deleting ${pvtu}"
        rm -f -- "$pvtu"

        echo "[delete] deleting partition files Output_p*_t${t_str}.vtu"
        find "$SIMDIR" -maxdepth 1 -type f -name "Output_p*_t${t_str}.vtu" -print -delete

    done < "$LIST"

    rm -f -- "$LIST"
    echo "[delete] Deleted list file: $LIST"
else
    echo "[delete] WARNING: list file missing or empty:"
    echo "$LIST"
fi

echo "Done."