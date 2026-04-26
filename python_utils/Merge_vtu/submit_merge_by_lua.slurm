#!/bin/bash
#SBATCH --job-name=Merge_by_lua
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

# Receive from sbatch --export (fallbacks if not provided)
SIMDIR="$ROOT"
SIMDIR="${SIMDIR%/}"   # Here is all the outputs (simulations)
T0="${T0:-0}"
T1="${T1:-0}"
SNAPSHOT_INTERVAL="${SNAPSHOT_INTERVAL:-15}"

PY_SCRIPT="/project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/Merge_vtu/merge_pvtu_to_vtu_by_lua.py"
LIST="${SIMDIR}/files_to_delete_T${T0}_T${T1}_J${SLURM_JOB_ID}.txt"

echo "[slurm] SIMDIR=$SIMDIR"
echo "[slurm] T0=$T0  T1=$T1  SNAPSHOT_INTERVAL=$SNAPSHOT_INTERVAL"


### STEP 1
time srun -n $_ntasks --ntasks-per-node=$_ntasks_per_node $HINT_FLAG -c $_nthreads --cpu-bind=verbose,cores  \
  `(which pvbatch)` --force-offscreen-rendering "$PY_SCRIPT" \
  --root "${SIMDIR}" \
  --list_to_delete "$LIST" \
  --t0 "${T0}" \
  --t1 "${T1}" \
  --snapshot_interval "${SNAPSHOT_INTERVAL}"


### STEP 2
#  delete fragments listed by python ----
if [[ -f "$LIST" ]]; then
  while IFS= read -r pvtu; do
    [[ -z "$pvtu" ]] && continue

    case "$pvtu" in
      "$SIMDIR"/Output_t*.pvtu) ;;
      *) echo "[delete] SKIP: $pvtu"; continue ;;
    esac

    t_str="$(basename "$pvtu")"
    t_str="${t_str#Output_t}"
    t_str="${t_str%.pvtu}"

    rm -f -- "$pvtu"
    find "$SIMDIR" -maxdepth 1 -type f -name "Output_p*_t${t_str}.vtu" -delete

  done < "$LIST"
  rm -f -- "$LIST"
else
  echo "[delete] WARNING: list file not found; nothing deleted"
fi