#!/bin/bash
EXEC=./executable
mkdir -p output logs

for case_dir in heterogeneity/*/; do
  echo "Starting simulation run '$case_dir' "
  case_name=$(basename "$case_dir")
  outdir="output/${case_name}"
  mkdir -p "$outdir"
  "$EXEC" "$case_dir" "$outdir/" > "logs/${case_name}.log" 2>&1
done
