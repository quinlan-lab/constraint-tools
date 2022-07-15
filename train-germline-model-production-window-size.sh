#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw

# https://devhints.io/bash#miscellaneous
# put option-fetching before "set -o nounset" so that we can detect flags without arguments
while [[ "$1" =~ ^- ]]; do 
  case $1 in
    --constraint-tools-directory ) shift; [[ ! $1 =~ ^- ]] && CONSTRAINT_TOOLS=$1;;
    --window-size ) shift; [[ ! $1 =~ ^- ]] && window_size=$1;;
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

progress_bars="disk" 
# progress_bars="stdout" 

work_directory="work-train-germline-model-production.windowSize-${window_size}"
work_directory_should_be_clean="true"
work="${CONSTRAINT_TOOLS_DATA}/${work_directory}" # path to directory to store intermediate work and logs
if [[ ${work_directory_should_be_clean} == "true" && -d ${work} ]]; then 
  error "the following work directory already exists:" ${work}
  exit 1 
else 
  mkdir --parents ${work}
fi 

trustworthy_noncoding_regions_train="${CONSTRAINT_TOOLS}/dist/trustworthy-noncoding-regions-germline-grch38-train.bed.gz"

model="${CONSTRAINT_TOOLS}/dist/model-germline-grch38.windowSize-${window_size}.json"

fetch_subset_of_trustworthy_noncoding_regions_train () { 
  less ${trustworthy_noncoding_regions_train} | head -10000
}

constraint-tools train-germline-model \
  --genome ${genome} \
  --build ${build} \
  --mutations ${mutations} \
  --number-chromosomes-min ${number_chromosomes_min} \
  --kmer-size ${kmer_size} \
  --model ${model} \
  --work ${work} \
  --progress-bars ${progress_bars} \
  --window-size ${window_size} \
  --trustworthy-noncoding-regions ${trustworthy_noncoding_regions_train}

  # TESTING:
  # --trustworthy-noncoding-regions <(fetch_subset_of_trustworthy_noncoding_regions_train | bgzip) \
  # --number-of-jobs 5 

