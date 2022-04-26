set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset

CONSTRAINT_TOOLS=$PWD

model="${CONSTRAINT_TOOLS}/dist/model-germline-grch38.json"
port="5000"

${CONSTRAINT_TOOLS}/constraint-tools dashboard-germline-model \
  --model ${model} \
  --port ${port}

  
