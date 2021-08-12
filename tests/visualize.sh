set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace

CONSTRAINT_TOOLS=$1

output="${CONSTRAINT_TOOLS}/tests" 
model="${output}/model.json"
port="5000"

${CONSTRAINT_TOOLS}/constraint-tools visualize \
  --model ${model} \
  --port ${port}

