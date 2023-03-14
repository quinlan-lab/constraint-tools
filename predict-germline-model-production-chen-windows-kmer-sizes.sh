#!/bin/bash
#SBATCH --time=40:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH --output=/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/genome-wide-predictions/predict-germline-model-production-chen-windows-kmer-sizes.log

set -o errexit
set -o pipefail
set -o nounset

source set-environment-variables.sh 

for kmer_size in "3" "5" "7"; do
  log_file="${CONSTRAINT_TOOLS_DATA}/genome-wide-predictions/predict-germline-grch38.chen-windows.kmerSize-${kmer_size}.log"
  job_name="predict-germline-model-production.chen-windows.kmerSize-${kmer_size}"
  sbatch \
    --wait \
    --output=${log_file} \
    --job-name=${job_name} \
  predict-germline-model-production-chen-windows-kmer-size.sh \
      --kmer-size ${kmer_size}
done
