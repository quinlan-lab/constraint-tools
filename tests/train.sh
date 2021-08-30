#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset
# set -o noclobber
# set -o xtrace

CONSTRAINT_TOOLS=$1 

<<<<<<< HEAD
kmer_size="3"
regions="${CONSTRAINT_TOOLS}/dist/neutral-regions-test.bed.gz" 
cell_type="germline"
model="${CONSTRAINT_TOOLS}/tests" # path to directory to store model in

if [ ${cell_type} == "somatic" ]; then
	mutations="${CONSTRAINT_TOOLS}/data/icgc/mutations.sorted.maf.gz"
	genome="${CONSTRAINT_TOOLS}/data/reference/grch37/genome.fa.gz"

elif [ ${cell_type} == "germline" ]; then
	mutations="${CONSTRAINT_TOOLS}/data/gnomad/v3/gnomad_v3_chr22.maf.gz"
	genome="${CONSTRAINT_TOOLS}/data/reference/grch38/hg38.analysisSet.fa.gz"

else 
	info "PLEASE SUPPLY \"germline\" OR \"somatic\" as input for the cell_type variable..."
fi
=======
mutations="/scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools/data/icgc/mutations.sorted.maf.gz"
genome="/scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools/data/reference/grch37/genome.fa.gz"
kmer_size="3"
regions="${CONSTRAINT_TOOLS}/tests/neutral-regions.bed.gz"
model="${CONSTRAINT_TOOLS}/tests" # path to directory to store model in
>>>>>>> upstream/main

${CONSTRAINT_TOOLS}/constraint-tools train \
  --genome ${genome} \
  --mutations ${mutations} \
  --kmer-size ${kmer_size} \
  --regions ${regions} \
<<<<<<< HEAD
  --cell-type ${cell_type} \
  --model ${model}
=======
  --model ${model}

>>>>>>> upstream/main
