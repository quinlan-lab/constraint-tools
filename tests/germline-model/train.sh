#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset
# set -o noclobber
# set -o xtrace

CONSTRAINT_TOOLS=$1 

genome="/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/reference/grch38/hg38.analysisSet.fa.gz"
mutations="/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/gnomad/v3/variants/gnomad_v3.sorted.tsv.gz"
number_chromosomes_min="120000"
kmer_size="3"
neutral_regions="${CONSTRAINT_TOOLS}/tests/germline-model/neutral-regions.bed.gz"
output="${CONSTRAINT_TOOLS}/tests/germline-model" 

${CONSTRAINT_TOOLS}/constraint-tools train-germline-model \
  --genome ${genome} \
  --mutations ${mutations} \
  --number-chromosomes-min ${number_chromosomes_min} \
  --kmer-size ${kmer_size} \
  --neutral-regions ${neutral_regions} \
  --output ${output}

