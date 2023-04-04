#!/bin/bash
#SBATCH --time=40:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH --output=/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/genome-wide-predictions/predict-germline-model-Nonly.production.chenWindows.kmerSizes.trainSets.log

set -o errexit
set -o pipefail
set -o nounset

source set-environment-variables.sh 

for kmer_size in "3" "5" "7"; do 
  # for train_set_label in "noncoding" "coding" "chenWindows"; do
  for train_set_label in "noncoding"; do
    log_file="${CONSTRAINT_TOOLS_DATA}/genome-wide-predictions/predict-germline-grch38-Nonly.chenWindows.kmerSize-${kmer_size}.trainSet-${train_set_label}.log"
    job_name="predict-germline-model-Nonly.production.chenWindows.kmerSize-${kmer_size}.trainSet-${train_set_label}"
    sbatch \
      --wait \
      --output=${log_file} \
      --job-name=${job_name} \
      predict-germline-model-Nonly.production.chenWindows.kmerSize.trainSet.sh \
        --kmer-size ${kmer_size} \
        --train-set-label ${train_set_label}
  done
done
