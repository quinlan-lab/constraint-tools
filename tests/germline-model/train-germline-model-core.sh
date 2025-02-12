#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH --output=dist/model-germline-grch38.log

all_trustworthy_noncoding_regions="true"
work_directory="work-train-germline-model"
work_directory_should_be_clean="true"

# https://devhints.io/bash#miscellaneous
# put option-fetching before "set -o nounset" so that we can detect flags without arguments
while [[ "$1" =~ ^- ]]; do 
  case $1 in
    --constraint-tools-directory ) shift; [[ ! $1 =~ ^- ]] && CONSTRAINT_TOOLS=$1;;
    --number-of-jobs ) shift; [[ ! $1 =~ ^- ]] && number_of_jobs=$1;;
    --number-of-trustworthy-noncoding-regions ) shift; [[ ! $1 =~ ^- ]] && number_of_trustworthy_noncoding_regions=$1 && all_trustworthy_noncoding_regions="false";;
    --work-directory ) shift; [[ ! $1 =~ ^- ]] && work_directory=$1;;
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
build="hg38" 
mutations="/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/gnomad/v3/variants/gnomad_v3.sorted.tsv.gz"
number_chromosomes_min="130000"
kmer_size="7"
work="${CONSTRAINT_TOOLS_DATA}/${work_directory}" # path to directory to store intermediate work and logs

progress_bars="disk" 
# progress_bars="stdout" 

if [[ ${work_directory_should_be_clean} == "true" && -d ${work} ]]; then 
  error "the following work directory already exists:" ${work}
  exit 1 
else 
  mkdir --parents ${work}
fi 

fetch_subset_of_trustworthy_noncoding_regions () { 
  local trustworthy_noncoding_regions="${CONSTRAINT_TOOLS}/dist/trustworthy-noncoding-regions-germline-grch38-train.bed.gz"
  less ${trustworthy_noncoding_regions} | head -${number_of_trustworthy_noncoding_regions}
}

train_on_subset_of_trustworthy_noncoding_regions () {
  # file to store model in:   
  local model="${CONSTRAINT_TOOLS}/tests/germline-model/model-germline-grch38-${number_of_trustworthy_noncoding_regions}.json" 

  # info "trustworthy noncoding regions (with their lengths):"
  # echo "$(fetch_subset_of_trustworthy_noncoding_regions | awk '{ print $0, $3-$2 }')"

  constraint-tools train-germline-model \
    --genome ${genome} \
    --build ${build} \
    --mutations ${mutations} \
    --number-chromosomes-min ${number_chromosomes_min} \
    --kmer-size ${kmer_size} \
    --trustworthy-noncoding-regions <(fetch_subset_of_trustworthy_noncoding_regions | bgzip) \
    --model ${model} \
    --work ${work} \
    --number-of-jobs ${number_of_jobs} \
    --progress-bars ${progress_bars}
}

train_on_all_trustworthy_noncoding_regions () {
  local model="${CONSTRAINT_TOOLS}/dist/model-germline-grch38.json" # file to store model in

  constraint-tools train-germline-model \
    --genome ${genome} \
    --build ${build} \
    --mutations ${mutations} \
    --number-chromosomes-min ${number_chromosomes_min} \
    --kmer-size ${kmer_size} \
    --model ${model} \
    --work ${work} \
    --progress-bars ${progress_bars}
}

if [[ ${all_trustworthy_noncoding_regions} == "true" ]]; then
  train_on_all_trustworthy_noncoding_regions
else
  train_on_subset_of_trustworthy_noncoding_regions
fi 
