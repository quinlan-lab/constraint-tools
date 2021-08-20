set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset

CONSTRAINT_TOOLS=$1

model="${CONSTRAINT_TOOLS}/tests/model.json" 
port="5000"

${CONSTRAINT_TOOLS}/constraint-tools visualize \
  --model ${model} \
  --port ${port}

  
