#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw

# https://devhints.io/bash#miscellaneous
# put option-fetching before "set -o nounset" so that we can detect flags without arguments
while [[ "$1" =~ ^- ]]; do 
  case $1 in
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

PATH=${CONSTRAINT_TOOLS}:$PATH 

model="${CONSTRAINT_TOOLS}/dist/model-germline-grch38.windowSize-${window_size}.json"
zscores="${CONSTRAINT_TOOLS}/dist/predict-germline-grch38.windowSize-${window_size}.bed.gz"
trustworthy_noncoding_regions="${CONSTRAINT_TOOLS}/dist/trustworthy-noncoding-regions-germline-grch38.bed.gz"
number_of_jobs="500"

progress_bars="disk" 
# progress_bars="stdout" 

work_directory="work-predict-germline-model-production.windowSize-${window_size}"
work_directory_should_be_clean="true"
work="${CONSTRAINT_TOOLS_DATA}/${work_directory}" # path to directory to store intermediate work and logs
if [[ ${work_directory_should_be_clean} == "true" && -d ${work} ]]; then 
  error "the following work directory already exists:" ${work}
  exit 1 
else 
  mkdir --parents ${work}
fi 

fetch_subset_of_trustworthy_noncoding_regions () { 
  less ${trustworthy_noncoding_regions} | head -10000 
}

constraint-tools predict-germline-model \
  --model ${model} \
  --zscores ${zscores} \
  --work ${work} \
  --progress-bars ${progress_bars} \
  --number-of-jobs ${number_of_jobs} \
  --trustworthy-noncoding-regions ${trustworthy_noncoding_regions}

  # TESTING:
  # --trustworthy-noncoding-regions <(fetch_subset_of_trustworthy_noncoding_regions | bgzip)

