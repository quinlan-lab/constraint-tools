#!/bin/bash

set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset

CONSTRAINT_TOOLS=$PWD

window_size="1001"
model="${CONSTRAINT_TOOLS}/dist/model-germline-grch38-exclude-test-promoters.windowSize-${window_size}.json"
port="5000"
trustworthy_noncoding_regions="${CONSTRAINT_TOOLS}/dist/trustworthy-noncoding-regions-germline-grch38.bed.gz"

"${CONSTRAINT_TOOLS}"/constraint-tools browse-germline-model \
  --model "${model}" \
  --port ${port} \
  --trustworthy-noncoding-regions "${trustworthy_noncoding_regions}"

  
