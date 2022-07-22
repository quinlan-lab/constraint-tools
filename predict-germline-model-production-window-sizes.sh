#!/bin/bash
#SBATCH --time=40:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH --output=dist/predict-germline-model-production-window-sizes.log

set -o errexit
set -o pipefail
set -o nounset

source set-environment-variables.sh 

for window_size in "101" "501" "1001" "2001" "5001"; do
  info "Computing predictions for model with window-size:" "${window_size}bp"
  log_file="dist/predict-germline-grch38.windowSize-${window_size}.log"
  job_name="predict-germline-model-production.windowSize-${window_size}"
  sbatch \
    --wait \
    --output=${log_file} \
    --job-name=${job_name} \
    predict-germline-model-production-window-size.sh \
      --window-size ${window_size}
done
info "Predictions complete for all window-sizes."
