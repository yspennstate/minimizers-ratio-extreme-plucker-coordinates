#!/bin/sh
# Runs every check and writes its output to outputs/<name>_out.txt.
set -e
PY=${PYTHON:-python3}
mkdir -p outputs
$PY -c "import sys, numpy, scipy; print('python', sys.version.split()[0], 'numpy', numpy.__version__, 'scipy', scipy.__version__)" > outputs/versions.txt
for s in ico lemma_check rank2_check; do
  $PY $s.py > outputs/${s}_out.txt
done
for s in ogr_check orbit_check col_check abs_check; do
  $PY $s.py outputs/${s}_out.txt > /dev/null
done
cat outputs/*_out.txt
