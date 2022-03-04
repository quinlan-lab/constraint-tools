#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH --output=dist/model-germline-grch38.log

all_neutral_regions="true"
work_subdirectory="work"
work_directory_should_be_clean="true"

# https://devhints.io/bash#miscellaneous
# put option-fetching before "set -o nounset" so that we can detect flags without arguments
while [[ "$1" =~ ^- ]]; do 
  case $1 in
    --constraint-tools-directory ) shift; [[ ! $1 =~ ^- ]] && CONSTRAINT_TOOLS=$1;;
    --number-of-neutral-regions ) shift; [[ ! $1 =~ ^- ]] && number_of_neutral_regions=$1 && all_neutral_regions="false";;
    --work-subdirectory ) shift; [[ ! $1 =~ ^- ]] && work_subdirectory=$1;;
    --work-directory-should-be-clean ) shift; [[ ! $1 =~ ^- ]] && work_directory_should_be_clean=$1;;
    *) error "$0: " "$1 is an invalid flag"; exit 1;;
  esac 
  shift
done

set -o errexit
set -o pipefail
set -o nounset
# set -o noclobber
# set -o xtrace

CONSTRAINT_TOOLS_DATA="/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools" 

PATH=${CONSTRAINT_TOOLS}:$PATH 
PATH=${CONSTRAINT_TOOLS}/utilities:$PATH 
PATH=${CONSTRAINT_TOOLS}/bin:$PATH 

genome="/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/reference/grch38/hg38.analysisSet.fa.gz"
mutations="/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/gnomad/v3/variants/gnomad_v3.sorted.tsv.gz"
number_chromosomes_min="130000"
kmer_size="3"
neutral_regions="${CONSTRAINT_TOOLS}/dist/neutral-regions-germline-grch38.bed.gz"
work="${CONSTRAINT_TOOLS_DATA}/${work_subdirectory}/train-germline-model" # path to directory to store intermediate work and logs

if [[ ${work_directory_should_be_clean} == "true" && -d ${work} ]]; then 
  error "the following work directory already exists:" ${work}
  exit 1 
else 
  mkdir --parents ${work}
fi 

fetch_subset_of_neutral_regions () { 
  less ${neutral_regions} | head -${number_of_neutral_regions}
}

train_on_subset_of_neutral_regions () {
  # file to store model in:   
  model="${CONSTRAINT_TOOLS}/tests/germline-model/model-germline-grch38-${number_of_neutral_regions}.json" 

  # 'stdout' or a directory to log progress to
  progress_bar="stdout" 

  # info "neutral regions (with their lengths):"
  # echo "$(fetch_subset_of_neutral_regions | awk '{ print $0, $3-$2 }')"

  constraint-tools train-germline-model \
    --genome ${genome} \
    --mutations ${mutations} \
    --number-chromosomes-min ${number_chromosomes_min} \
    --kmer-size ${kmer_size} \
    --neutral-regions <(fetch_subset_of_neutral_regions | bgzip) \
    --model ${model} \
    --work ${work} \
    --number-of-jobs "2" \
    --progress-bar ${progress_bar}
}

train_on_all_neutral_regions () {
  model="${CONSTRAINT_TOOLS}/dist/model-germline-grch38.json" # file to store model in

  info "training on all neutral regions..."

  constraint-tools train-germline-model \
    --genome ${genome} \
    --mutations ${mutations} \
    --number-chromosomes-min ${number_chromosomes_min} \
    --kmer-size ${kmer_size} \
    --neutral-regions ${neutral_regions} \
    --model ${model} \
    --work ${work} \
    --progress-bar ${progress_bar}
}

if [[ ${all_neutral_regions} == "true" ]]; then
  train_on_all_neutral_regions
else
  train_on_subset_of_neutral_regions
fi 
