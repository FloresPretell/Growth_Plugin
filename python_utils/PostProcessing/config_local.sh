#!/usr/bin/env bash

DATASET_ROOT="/Users/flore0a/Downloads/Fusion_pipeline/PostProcessing/Cases"
OUTPUT_ROOT="/Users/flore0a/Downloads/Fusion_pipeline/PostProcessing/Results"
MAX_CASES=""

# Optional local Conda environment for ParaView + Python tools.
# Example:
#   conda activate paraview_env
#   which python
#   which pvpython
#   which pvbatch
#
# If you prefer absolute paths, set the interpreter variables below directly.
LOCAL_CONDA_ENV="paraview_env"

LOCAL_PYTHON="python"
PNG_EXPORT_PYTHON="${LOCAL_PYTHON}"
PNG_EXPORT_MODE="local_safe"

# Local quick-test defaults. Clear these values for a full local_safe export.
# Stage 2 and stage 3 consume u_t PNGs, so limiting local testing to u_t keeps
# the downstream chain testable without matching the full Shaheen workload.
LOCAL_FIELDS_TO_EXPORT="u_t"
LOCAL_MAX_TIMESTEPS="1"

MATRIX_PYTHON="${LOCAL_PYTHON}"
# Local SWC/morphology stages need packages that are not present in paraview_env.
SWC_PYTHON="/opt/anaconda3/envs/swc_env/bin/python"
MORPHO_PYTHON="/opt/anaconda3/envs/swc_env/bin/python"

SLURM_ACCOUNT=""
SLURM_PARTITION=""
SLURM_MAIL=""
PVPYTHON_MODULE=""
PVPYTHON_SOURCE=""
