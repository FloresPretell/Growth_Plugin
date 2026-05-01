#!/bin/bash
set -euo pipefail

ENV_NAME="swc_env"
PYTHON_VERSION="3.10"

echo "=== Loading conda ==="

# Adjust this block if your Shaheen conda path is different
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
elif command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
else
    echo "ERROR: conda not found. Load your conda installation first."
    exit 1
fi

echo "=== Checking whether environment ${ENV_NAME} already exists ==="
if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
    echo "Environment ${ENV_NAME} already exists."
else
    echo "Creating environment ${ENV_NAME} with Python ${PYTHON_VERSION}"
    conda create -y -n "${ENV_NAME}" python="${PYTHON_VERSION}"
fi

echo "=== Activating ${ENV_NAME} ==="
conda activate "${ENV_NAME}"

echo "=== Installing packages ==="
conda install -y -c conda-forge \
    numpy \
    scipy \
    pillow \
    scikit-image \
    matplotlib \
    jupyter \
    ipykernel

echo "=== Verifying installation ==="
python - <<'PY'
import sys
import numpy
import scipy
import skimage
from PIL import Image
import matplotlib

print("Python executable:", sys.executable)
print("numpy:", numpy.__version__)
print("scipy:", scipy.__version__)
print("skimage:", skimage.__version__)
print("Pillow OK")
print("matplotlib:", matplotlib.__version__)
print("Environment installation successful.")
PY

echo "=== Done. Environment ${ENV_NAME} is ready. ==="