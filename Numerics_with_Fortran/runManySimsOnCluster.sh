#!/bin/bash
#SBATCH --job-name=sim_ensemble
#SBATCH --partition=geo4dgpu
#SBATCH --gres=gpu:1
#SBATCH --mem=12G
#SBATCH --time=08:00:00
#SBATCH --output=logs/slurm_master.out

module load nvhpc-hpcx-cuda12/25.3   # module loading for nvhpc
module load cuda/12.8.1

EXEC=./executable
mkdir -p output logs

for case_dir in heterogeneity/*/; do
  case_name=$(basename "$case_dir")
  outdir="output/${case_name}"
  mkdir -p "$outdir"

  echo "Starting simulation run '$case_name'"
  "$EXEC" "$case_dir" "$outdir/" > >(tee "logs/${case_name}.log") 2>&1

  if [ $? -ne 0 ]; then
    echo "  FAILED — see logs/${case_name}.log"
  else
    echo "  OK"
  fi
done
