#!/bin/bash
#SBATCH --time=20:00:00
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
    --directory ) shift; [[ ! $1 =~ ^- ]] && directory=$1;;
    --window-index ) shift; [[ ! $1 =~ ^- ]] && window_index=$1;;
    *) error "$0: " "$1 is an invalid flag"; exit 1;;
  esac 
  shift
done

set -o errexit
set -o pipefail
set -o nounset
# set -o noclobber
# set -o xtrace

PATH="${CONSTRAINT_TOOLS}/experiments/germline-model/chen-et-al-2022:$PATH" 

info "Inferring using nucleotide transformer..."
info "\tdirectory:" ${directory}
info "\twindow index:" ${window_index}

nucleotide-transformer.2.inference "${directory}" "${window_index}"

