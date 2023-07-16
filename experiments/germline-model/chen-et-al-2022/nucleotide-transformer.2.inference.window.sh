#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-rw
#SBATCH --mem=50g # sacct -o reqmem,maxrss,averss,elapsed -j JOBID

# $ myallocation

# Redwood general cluster hardware overview: 
# https://www.chpc.utah.edu/documentation/guides/redwood.php#apexHardware

# previously used as much as: 
# --mem=192g

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

