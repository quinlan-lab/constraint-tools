set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset

CONSTRAINT_TOOLS=$PWD

window_size="1001"
model="${CONSTRAINT_TOOLS}/dist/model-germline-grch38-exclude-test-promoters.windowSize-${window_size}.json"
port="5000"

${CONSTRAINT_TOOLS}/constraint-tools browser-germline-model \
  --model ${model} \
  --port ${port}

  
