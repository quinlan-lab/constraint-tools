#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH --mem=20g # sacct -o reqmem,maxrss,averss,elapsed -j JOBID

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
genome_wide_predictions_directory="${CONSTRAINT_TOOLS_DATA}/genome-wide-predictions"
zscores="${genome_wide_predictions_directory}/predict-germline-grch38.windowSize-${window_size}.bed.gz"

progress_bars="disk" 
# progress_bars="stdout" 

work_directory="work-predict-germline-model-production.windowSize-${window_size}"
work_directory_should_be_clean="true"
work="${genome_wide_predictions_directory}/${work_directory}" # path to directory to store intermediate work and logs
if [[ ${work_directory_should_be_clean} == "true" && -d ${work} ]]; then 
  error "the following work directory already exists:" ${work}
  exit 1 
else 
  mkdir --parents ${work}
fi 

trustworthy_noncoding_windows="${work}/trustworthy-noncoding-windows.bed" 
create-trustworthy-noncoding-windows \
  --model ${model} \
  --work ${work} \
  --output ${trustworthy_noncoding_windows}

# TESTING: 
# head -485 ${trustworthy_noncoding_windows} > trustworthy_noncoding_windows_small.bed

constraint-tools predict-germline-model \
  --model ${model} \
  --zscores ${zscores} \
  --work ${work} \
  --progress-bars ${progress_bars} \
  --number-of-jobs 500 \
  --windows ${trustworthy_noncoding_windows}

# TESTING: 
  # --number-of-jobs 5 \
  # --windows trustworthy_noncoding_windows_small.bed


