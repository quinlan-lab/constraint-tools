#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH --output=experiments/germline-model/cpg-islands/cpg-islands.log

jupyter nbconvert \
  --execute \
  --to notebook \
  --inplace \
  --ExecutePreprocessor.timeout -1 \
  --NbConvertApp.log_level DEBUG \
  experiments/germline-model/cpg-islands/cpg-islands.ipynb
