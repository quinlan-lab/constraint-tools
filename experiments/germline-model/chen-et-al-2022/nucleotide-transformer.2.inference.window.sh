#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-rw
#SBATCH --mem=50g # sacct -o reqmem,maxrss,averss,elapsed -j JOBID

######################## [START] Notes on memory requirements #################################

# "The model was trained with 8 A100 80GB" : 
# https://huggingface.co/InstaDeepAI/nucleotide-transformer-500m-human-ref
# Back of the envelope calculations of memory requirement of transformer models: 
# https://twitter.com/MishaLaskin/status/1546994229674647553?s=20&t=0gkdvE1j_363D3xvTT1d4A
# "With commonly available current hardware and model sizes, 
# [computational and memory requirements that are quadratic with the input sequence length] 
# typically limits the input sequence to roughly 512 tokens" : 
# https://ai.googleblog.com/2021/03/constructing-transformers-for-longer.html

# $ myallocation

# could use as much as: 
# --mem=192g
# on redwood: 
# https://www.chpc.utah.edu/documentation/guides/redwood.php#apexHardware

######################## [END] Notes on memory requirements #################################

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

