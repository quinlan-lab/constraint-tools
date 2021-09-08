#!/bin/bash

#SBATCH --time=3:00:00
#SBATCH --nodes=4
#SBATCH --ntasks=16
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH -o logs/train-model.out

set -o errexit
set -o pipefail
set -o nounset
# set -o noclobber
# set -o xtrace

CONSTRAINT_TOOLS=$1 

kmer_size="3"
regions="${CONSTRAINT_TOOLS}/dist/neutral-regions-test.bed.gz" 
cell_type="germline"
model="${CONSTRAINT_TOOLS}/tests" # path to directory to store model in

if [ ${cell_type} == "somatic" ]; then
	mutations="${CONSTRAINT_TOOLS}/data/icgc/mutations.sorted.maf.gz"
	genome="${CONSTRAINT_TOOLS}/data/reference/grch37/genome.fa.gz"

elif [ ${cell_type} == "germline" ]; then
	mutations="${CONSTRAINT_TOOLS}/data/gnomad/v3/gnomad_v3_variants.sorted.bed.gz"
	genome="${CONSTRAINT_TOOLS}/data/reference/grch38/hg38.analysisSet.fa.gz"

else 
	info "PLEASE SUPPLY \"germline\" OR \"somatic\" as input for the cell_type variable..."
fi

kmer_size="3"
regions="${CONSTRAINT_TOOLS}/tests/neutral-regions.bed.gz"
model="${CONSTRAINT_TOOLS}/tests" # path to directory to store model in

${CONSTRAINT_TOOLS}/constraint-tools train \
  --genome ${genome} \
  --mutations ${mutations} \
  --kmer-size ${kmer_size} \
  --regions ${regions} \
  --cell-type ${cell_type} \
  --model ${model}
