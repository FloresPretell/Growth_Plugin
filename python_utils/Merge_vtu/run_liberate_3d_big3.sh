#!/bin/bash
source /scratch/flore0a/Modules.sh mesa >/dev/null 2>&1
PVB=/sw/vis2/shaheen3.2026.5/paraview-buildMesa-6.1.0/install/bin/pvbatch
MV=/project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/Merge_vtu/merge_and_purge_3d.py
BASE=/scratch/flore0a/3D_evaluator
inodes(){ lfs quota -u flore0a /scratch 2>/dev/null | awk '/\/scratch/{print $6}'; }
echo "START inodes=$(inodes)  $(date)"
DIRS=$(find $BASE/realeik2_20260706_12749499 $BASE/eikonaldir_20260706_12749356 $BASE/eikref_20260706_12749511 -type d -name growth 2>/dev/null | sort)
i=0
for D in $DIRS; do
  ls $D/*_t*.pvtu >/dev/null 2>&1 || continue
  i=$((i+1))
  echo "### [$i] $D  ($(date))"
  $PVB $MV --root "$D" 2>&1 | grep -v -iE 'warning|vtkmodules|generic'
  echo "   inodes now: $(inodes)"
done
echo "=== ALL DONE $(date)  inodes=$(inodes) ==="
