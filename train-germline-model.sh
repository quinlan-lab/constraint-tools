#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH --output=dist/train-germline-model.log

set -o errexit
set -o pipefail
set -o nounset
# set -o noclobber
# set -o xtrace

CONSTRAINT_TOOLS=$1 
CONSTRAINT_TOOLS_DATA="/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools" 

PATH=${CONSTRAINT_TOOLS}:$PATH 
PATH=${CONSTRAINT_TOOLS}/utilities:$PATH 
PATH=${CONSTRAINT_TOOLS}/bin:$PATH 

genome="/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/reference/grch38/hg38.analysisSet.fa.gz"
mutations="/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/gnomad/v3/variants/gnomad_v3.sorted.tsv.gz"
number_chromosomes_min="130000"
kmer_size="3"
neutral_regions="${CONSTRAINT_TOOLS}/dist/neutral-regions-germline-grch38.bed.gz"
model="${CONSTRAINT_TOOLS}/dist/model-germline-grch38.json" # file to store model in
work="${CONSTRAINT_TOOLS_DATA}/work/train-germline-model" # path to directory to store intermediate work and logs
mkdir --parents ${work}

fetch_subset_of_neutral_regions () { 
  less ${neutral_regions} | head -100 
}

train_on_subset_of_neutral_regions () {
  info "neutral regions (with their lengths):"
  echo "$(fetch_subset_of_neutral_regions | awk '{ print $0, $3-$2 }')"

  constraint-tools train-germline-model \
    --genome ${genome} \
    --mutations ${mutations} \
    --number-chromosomes-min ${number_chromosomes_min} \
    --kmer-size ${kmer_size} \
    --neutral-regions <(fetch_subset_of_neutral_regions | bgzip) \
    --model ${model} \
    --work ${work}
}

train_on_all_neutral_regions () {
  info "training on all neutral regions..."

  constraint-tools train-germline-model \
    --genome ${genome} \
    --mutations ${mutations} \
    --number-chromosomes-min ${number_chromosomes_min} \
    --kmer-size ${kmer_size} \
    --neutral-regions ${neutral_regions} \
    --model ${model} \
    --work ${work}
}

train_on_subset_of_neutral_regions
# train_on_all_neutral_regions