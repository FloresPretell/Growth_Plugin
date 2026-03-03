#!/bin/bash
#SBATCH -N 1                         # 8 nodos
#SBATCH -J small_test                   # nombre del job
#SBATCH --partition=workq
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
#SBATCH --mail-user=elsanicole.florespretell@kaust.edu.sa
#SBATCH --mail-type=ALL
#SBATCH -A k10070
#SBATCH -t 23:40:00

mkdir -p logs

#OpenMP settings:
export OMP_NUM_THREADS=1

directory="/scratch/flore0a/Dataset_test1/small_test"
mkdir -p "$directory"
cd "$directory"

# --- Start timer ---
start_ts=$(date +%s)

# ==== 1) Simulation ====
srun --hint=nomultithread \
     --ntasks=96 \
     /project/k10070/Nicole/UG4/ug4/bin/ugshell \
       -ex /project/k10070/Nicole/UG4/ug4/apps/Neuronal_Process/Scripts_lua/Model_only_final_steps.lua \
       -grid gradiente4.ugx \
       -inhibition \
       -inside \
       -gradiente1 \
       -gradiente2 \
       -deformation \
       -reaction1 -reaction2 -reaction3 -reaction4

# --- End timer ---
end_ts=$(date +%s)
elapsed=$(( end_ts - start_ts ))
# --- Save timing ---
echo "time=${elapsed}s" >> timing_test2.txt

# ==== 2) Submit merge job with dependency ====
PACK_SLURM="/project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/Merge_vtu/submit_merge_one_file_by_simulation.slurm"
CASE_ROOT="${directory}/Simulation2D"   # aquí está merged/

# Solo si existen VTUs mergeados
if ls "${CASE_ROOT}/merged/Output_t"*.vtu >/dev/null 2>&1; then
  sbatch --dependency=afterok:${SLURM_JOB_ID} \
         --export=ALL,ROOT="${CASE_ROOT}" \
         "$PACK_SLURM"
else
  echo "[pack] No merged VTU files found in ${CASE_ROOT}/merged, skipping pack."
fi