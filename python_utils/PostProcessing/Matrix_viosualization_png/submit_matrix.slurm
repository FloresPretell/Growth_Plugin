#!/bin/bash
#SBATCH --job-name=Export_swc
#SBATCH --partition=workq
#SBATCH --mail-user=elsanicole.florespretell@kaust.edu.sa
#SBATCH --mail-type=ALL
#SBATCH -A k10070
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --array=0-1
#SBATCH --mem=32G
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail
PYTHON_BIN="/scratch/flore0a/iops/miniconda3/envs/swc_env/bin/python"

mkdir -p logs

SIMDIR=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" dirs.txt)

if [ -z "${SIMDIR}" ]; then
    echo "ERROR: No directory found in dirs.txt for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "Task ${SLURM_ARRAY_TASK_ID} exporting SWC from ${SIMDIR}"

"${PYTHON_BIN}" create_mosaic_gif_matrix.py --base_dir "${SIMDIR}" --output "${SIMDIR}/mosaic_matrix.gif" --tile_size 140 100 --duration 250 --keep_last_frame --png_pattern "*.png" --png_subdir "Simulation2D/merged/png/u_t"



