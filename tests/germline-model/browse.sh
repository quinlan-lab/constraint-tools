#!/bin/bash

set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset

CONSTRAINT_TOOLS=$PWD

window_size="1000"
model="${CONSTRAINT_TOOLS}/dist/model-germline-grch38.windowSize-${window_size}.json"
port="5001" # https://stackoverflow.com/a/69829313
trustworthy_noncoding_regions="${CONSTRAINT_TOOLS}/dist/trustworthy-noncoding-regions-germline-grch38.bed.gz"

"${CONSTRAINT_TOOLS}"/constraint-tools browse-germline-model \
  --port ${port} \
  --model "${model}" \
  --trustworthy-noncoding-regions "${trustworthy_noncoding_regions}"

  
