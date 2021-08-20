#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset
# set -o noclobber
# set -o xtrace

CONSTRAINT_TOOLS=$1 

mutations="/scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools/data/icgc/mutations.sorted.maf.gz"
genome="/scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools/data/reference/grch37/genome.fa.gz"
kmer_size="3"
regions="${CONSTRAINT_TOOLS}/tests/neutral-regions.bed.gz"
model="${CONSTRAINT_TOOLS}/tests" # path to directory to store model in

${CONSTRAINT_TOOLS}/constraint-tools train \
  --genome ${genome} \
  --mutations ${mutations} \
  --kmer-size ${kmer_size} \
  --regions ${regions} \
  --model ${model}

