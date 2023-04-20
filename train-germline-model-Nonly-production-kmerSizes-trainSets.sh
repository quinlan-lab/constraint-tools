#!/bin/bash
#SBATCH --time=40:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH --output=dist/train-germline-model-Nonly-production-kmerSizes-trainSets.log

set -o errexit 
set -o nounset

CONSTRAINT_TOOLS="${PWD}"
CONSTRAINT_TOOLS_DATA="/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools" 

trustworthy_noncoding_regions="${CONSTRAINT_TOOLS}/dist/trustworthy-noncoding-regions-germline-grch38.bed.gz"
merged_exons="${CONSTRAINT_TOOLS_DATA}/genes/grch38/exons.merged.bed.gz"
chen_windows="${CONSTRAINT_TOOLS_DATA}/chen-et-al-2022/chen-windows.bed.gz"

# for kmer_size in "3" "5" "7"; do 
for kmer_size in "1"; do 
  for trainSet_trainSetLabel in \
      "${trustworthy_noncoding_regions},noncoding"; do
      # "${merged_exons},coding" \
      # "${chen_windows},chenWindows"; do
    IFS=, read train_set train_set_label <<< ${trainSet_trainSetLabel}
    log_file="dist/model-germline-grch38-Nonly.kmerSize-${kmer_size}.trainSet-${train_set_label}.log"
    job_name="train-germline-model-Nonly-production.kmerSize-${kmer_size}.trainSet-${train_set_label}"
    sbatch \
      --wait \
      --output=${log_file} \
      --job-name=${job_name} \
      train-germline-model-Nonly-production-kmerSize-trainSet.sh \
        --constraint-tools-directory ${CONSTRAINT_TOOLS} \
        --constraint-tools-data-directory ${CONSTRAINT_TOOLS_DATA} \
        --kmer-size ${kmer_size} \
        --train-set ${train_set} \
        --train-set-label ${train_set_label}
  done
done
