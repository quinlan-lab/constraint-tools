#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --nodes=2
#SBATCH --ntasks=16
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw

set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace

CONSTRAINT_TOOLS=$1 

mutations="/scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools/data/icgc/mutations.sorted.maf.gz"
genome="/scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools/data/reference/grch37/genome.fa.gz"
kmer_size="5"
output="${CONSTRAINT_TOOLS}/tests" 

${CONSTRAINT_TOOLS}/constraint-tools train \
  --genome ${genome} \
  --mutations ${mutations} \
  --kmer-size ${kmer_size} \
  --output ${output}

