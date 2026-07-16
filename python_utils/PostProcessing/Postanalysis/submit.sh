#!/bin/bash

set -euo pipefail
directorio=/scratch/flore0a/Dataset_test2/Cases

# 1) Generar lista de directorios
python3 ./generate_folder_list.py --path "$directorio"

# 2) Contar cuántas tareas
N=$(wc -l < dirs.txt)
if [[ "$N" -le 0 ]]; then
  echo "ERROR: dirs.txt está vacío"
  exit 1
fi
LAST=$((N-1))

echo "Found $N simulations -> submitting array 0-$LAST"

# 3)Sent SLURM, sobrescribiendo el --array del slurm script
sbatch --array=0-"$LAST" ./Submite_export.slurm