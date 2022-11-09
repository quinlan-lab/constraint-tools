#!/bin/bash
#SBATCH --time=40:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH --output=dist/train-germline-model-production-window-sizes.log

for window_size in "100" "500" "1000" "2000" "5000"; do
  log_file="dist/model-germline-grch38.windowSize-${window_size}.log"
  job_name="train-germline-model-production.windowSize-${window_size}"
  sbatch \
    --wait \
    --output=${log_file} \
    --job-name=${job_name} \
    train-germline-model-production-window-size.sh \
      --constraint-tools-directory $PWD \
      --window-size ${window_size}
done
