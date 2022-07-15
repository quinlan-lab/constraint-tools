#!/bin/bash
#SBATCH --time=40:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH --output=dist/train-germline-model-production-window-sizes.log

for window_size in "101" "501" "1001" "2001" "5001"; do
  log_file="dist/model-germline-grch38-exclude-test-promoters.windowSize-${window_size}.log"
  job_name="train-germline-model-production.windowSize-${window_size}"
  sbatch \
    --wait \
    --output=${log_file} \
    --job-name=${job_name} \
    train-germline-model.exclude-test-promoters.sh \
      --constraint-tools-directory $PWD \
      --window-size ${window_size}
done
