#!/usr/bin/env bash

# set -o xtrace 

model="${CONSTRAINT_TOOLS}/$(read-config predictGermlineModel model)"
trustworthy_noncoding_regions="${CONSTRAINT_TOOLS}/$(read-config predictGermlineModel trustworthyNoncodingRegions)"

# https://devhints.io/bash#miscellaneous
# put option-fetching before "set -o nounset" so that we can detect flags without arguments
while [[ "$1" =~ ^- ]]; do 
  case $1 in
    --model ) shift; [[ ! $1 =~ ^- ]] && model=$1;;
    --trustworthy-noncoding-regions ) shift; [[ ! $1 =~ ^- ]] && trustworthy_noncoding_regions=$1;;
    *) error "$0: $1 is an invalid flag"; exit 1;;
  esac 
  shift
done

set -o errexit
set -o pipefail
set -o noclobber
set -o nounset

info "Using the model specified at:" ${model}
info "Using the trustworthy noncoding regions specified at:" ${trustworthy_noncoding_regions}

# TODO: 
# 1. tile the genome into window_size tiles 
# 2. find tiles that lie within trustworthy_noncoding_regions 
# 3. divide tiles into 500 batches
# 4. MAP | for each batch, launch a python program that loops over tiles in batch and dumps scores for all tiles to disk 
#    predict-constraint/germline-model/expected_observed_counts.py > compute_zscores     
# 5. REDUCE | concatenate all score files 
# 6. uupdate README to describe usage of this script
