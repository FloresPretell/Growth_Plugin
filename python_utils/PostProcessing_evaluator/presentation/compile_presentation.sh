#!/bin/bash
# compile_presentation.sh — build presentation.pdf from presentation.tex.
#
# pdflatex is provided by the Shaheen system TeX Live 2022 installation; the
# binary is not in the default PATH, so point at it directly. Two passes
# resolve table-of-contents / hyperref cross-references. Same pattern as
# Promesh_evaluator/presentation/compile_presentation.sh.
set -euo pipefail
cd "$(dirname "$0")"

PDFLATEX="/sw/rl9c/texlive/2022/rl9_binary/install-tl-20221109/bin/x86_64-linux/pdflatex"
[[ -x "${PDFLATEX}" ]] || { echo "ERROR: pdflatex not found at ${PDFLATEX}"; exit 2; }

"${PDFLATEX}" -interaction=nonstopmode -halt-on-error presentation.tex
"${PDFLATEX}" -interaction=nonstopmode -halt-on-error presentation.tex

echo
echo "Build complete."
ls -lh presentation.pdf
