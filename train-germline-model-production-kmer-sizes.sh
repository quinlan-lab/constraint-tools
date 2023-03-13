#!/bin/bash
#SBATCH --time=40:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH --output=dist/train-germline-model-production-kmer-sizes.log

for kmer_size in "3" "5" "7"; do
  log_file="dist/model-germline-grch38.kmerSize-${kmer_size}.log"
  job_name="train-germline-model-production.kmerSize-${kmer_size}"
  sbatch \
    --wait \
    --output=${log_file} \
    --job-name=${job_name} \
    train-germline-model-production-kmer-size.sh \
      --constraint-tools-directory $PWD \
      --kmer-size ${kmer_size}
done
