#!/bin/bash
#SBATCH --time=40:00:00
#SBATCH --account=redwood-gpu
#SBATCH --partition=redwood-gpu
#SBATCH --output=/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/genome-wide-predictions/predict-germline-model-production-window-sizes.log

# previously used: 
# --account=quinlan-rw
# --partition=quinlan-shared-rw

set -o errexit
set -o pipefail
set -o nounset

source set-environment-variables.sh 

for window_size in "100" "500" "1000" "2000" "5000" "10000"; do
  info "Computing predictions for model with window-size:" "${window_size}bp"
  log_file="${CONSTRAINT_TOOLS_DATA}/genome-wide-predictions/predict-germline-grch38.windowSize-${window_size}.log"
  job_name="predict-germline-model-production.windowSize-${window_size}"
  sbatch \
    --wait \
    --output=${log_file} \
    --job-name=${job_name} \
    predict-germline-model-production-window-size.sh \
      --window-size ${window_size}
done
info "Predictions complete for all window-sizes."
