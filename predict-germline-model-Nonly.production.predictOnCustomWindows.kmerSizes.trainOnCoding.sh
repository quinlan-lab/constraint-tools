#!/bin/bash
#SBATCH --time=40:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH --output=/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/genome-wide-predictions/predict-germline-model-Nonly.production.predictOnCustomWindows.kmerSizes.trainOnCoding.log

set -o errexit
set -o pipefail
set -o nounset

source set-environment-variables.sh 

windows_to_predict_on=${PUT_PATH_TO_YOUR_WINDOWS_HERE} 

for kmer_size in "3" "5" "7"; do 
  log_file="${CONSTRAINT_TOOLS_DATA}/genome-wide-predictions/predict-germline-grch38-Nonly.predictOnCustomWindows.kmerSize-${kmer_size}.trainOnCoding.log"
  job_name="predict-germline-model-Nonly.production.predictOnCustomWindows.kmerSize-${kmer_size}.trainOnCoding"
  sbatch \
    --wait \
    --output=${log_file} \
    --job-name=${job_name} \
    predict-germline-model-Nonly.production.predictOnCustomWindows.kmerSize.trainOnCoding.sh \
      --kmer-size ${kmer_size} \
      --windows-to-predict-on ${windows_to_predict_on}
done
