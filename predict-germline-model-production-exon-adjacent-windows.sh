#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

source set-environment-variables.sh 

log_file="${CONSTRAINT_TOOLS_DATA}/genome-wide-predictions/predict-germline-grch38.exon-adjacent-windows.log"
job_name="predict-germline-model-production.exon-adjacent-windows"
sbatch \
  --output=${log_file} \
  --job-name=${job_name} \
  predict-germline-model-production-exon-adjacent-windows-core.sh
