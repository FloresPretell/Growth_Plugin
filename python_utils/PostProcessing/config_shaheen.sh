#!/usr/bin/env bash

DATASET_ROOT="/scratch/flore0a/Dataset_test4/Cases"
OUTPUT_ROOT="/scratch/flore0a/Dataset_test4_Results"
MAX_CASES=""

LOCAL_PYTHON="python3"
PVPYTHON_CMD="pvbatch"
PVPYTHON_ARGS="--force-offscreen-rendering"
PNG_EXPORT_MODE="render"
PVPYTHON_SOURCE="/scratch/flore0a/Modules.sh"

# === Use conda environment for all non-ParaView tasks ===
CONDA_PYTHON="/scratch/flore0a/iops/miniconda3/envs/swc_env/bin/python"

MATRIX_PYTHON="${CONDA_PYTHON}"      # ✅ Use conda
SWC_PYTHON="${CONDA_PYTHON}"         # ✅ Use conda  
MORPHO_PYTHON="${CONDA_PYTHON}"      # ✅ Use conda

SLURM_ACCOUNT="k10070"
SLURM_PARTITION="workq"
SLURM_MAIL="elsanicole.florespretell@kaust.edu.sa"