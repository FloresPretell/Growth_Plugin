#!/bin/bash
source /scratch/flore0a/Modules.sh mesa >/dev/null 2>&1
PVB=/sw/vis2/shaheen3.2026.5/paraview-buildMesa-6.1.0/install/bin/pvbatch
MV=/project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/Merge_vtu/merge_and_purge_3d.py
LIST=/tmp/claude-194965/-scratch-flore0a/6bb62a79-db00-4d29-b3cd-d11bf40ec0c6/scratchpad/turing_prep_dirs.txt
inodes(){ lfs quota -u flore0a /scratch 2>/dev/null | awk '/\/scratch/{print $6}'; }
n=$(wc -l < "$LIST"); i=0
echo "START turing+prep inodes=$(inodes)  $(date)  ($n dirs)"
while IFS= read -r D; do
  [ -d "$D" ] || continue
  ls "$D"/*_t*.pvtu >/dev/null 2>&1 || continue
  i=$((i+1)); echo "### [$i/$n] $D"
  $PVB $MV --root "$D" 2>&1 | grep -v -iE 'warning|vtkmodules|generic' | tail -2
  echo "   inodes now: $(inodes)"
done < "$LIST"
echo "=== TURING+PREP DONE $(date)  inodes=$(inodes) ==="
