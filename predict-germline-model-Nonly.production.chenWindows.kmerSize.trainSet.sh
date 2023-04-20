#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --account=redwood-gpu
#SBATCH --partition=redwood-gpu
#SBATCH --mem=20g # sacct -o reqmem,maxrss,averss,elapsed -j JOBID

# previously used: 
# --account=quinlan-rw
# --partition=quinlan-shared-rw

# https://devhints.io/bash#miscellaneous
# put option-fetching before "set -o nounset" so that we can detect flags without arguments
while [[ "$1" =~ ^- ]]; do 
  case $1 in
    --kmer-size ) shift; [[ ! $1 =~ ^- ]] && kmer_size=$1;;
    --train-set-label ) shift; [[ ! $1 =~ ^- ]] && train_set_label=$1;;
    --noisy-label ) shift; [[ ! $1 =~ ^- ]] && noisy_label=$1;;
    *) error "$0: " "$1 is an invalid flag"; exit 1;;
  esac 
  shift
done

set -o errexit
set -o pipefail
set -o nounset
# set -o noclobber
# set -o xtrace

info "Computing Nbars on Chen windows using:"
info "\tkmer size:" ${kmer_size}bp
info "\ttraining set:" ${train_set_label}
info "\tnoisy label:" ${noisy_label}

PATH=${CONSTRAINT_TOOLS}:$PATH 

model="${CONSTRAINT_TOOLS}/dist/model-germline-grch38-Nonly.kmerSize-${kmer_size}.trainSet-${train_set_label}${noisy_label}.json"

genome_wide_predictions_directory="${CONSTRAINT_TOOLS_DATA}/genome-wide-predictions"
zscores="${genome_wide_predictions_directory}/predict-germline-grch38-Nonly.chenWindows.kmerSize-${kmer_size}.trainSet-${train_set_label}${noisy_label}.bed.gz"

progress_bars="disk" 
# progress_bars="stdout" 

work_directory="work-predict-germline-model-Nonly.production.chenWindows.kmerSize-${kmer_size}.trainSet-${train_set_label}${noisy_label}"
work_directory_should_be_clean="true"
work="${genome_wide_predictions_directory}/${work_directory}" # path to directory to store intermediate work and logs
if [[ ${work_directory_should_be_clean} == "true" && -d ${work} ]]; then 
  error "the following work directory already exists:" ${work}
  exit 1 
else 
  mkdir --parents ${work}
fi 

chen_windows_with_scores="${CONSTRAINT_TOOLS_DATA}/chen-et-al-2022/Supplementary_Datasets/Supplementary_Data_2.bed"
chen_windows="${work}/chen-windows.bed" 
cut -f1-3 ${chen_windows_with_scores} | tail -n +2 > ${chen_windows}
info "Windows:" ${chen_windows}

constraint-tools predict-germline-model-Nonly \
  --model ${model} \
  --zscores ${zscores} \
  --work ${work} \
  --progress-bars ${progress_bars} \
  --number-of-jobs 500 \
  --windows ${chen_windows}


