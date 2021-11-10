set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset

CONSTRAINT_TOOLS=$1

# model="${CONSTRAINT_TOOLS}/tests/model.json" 
model="/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/model.ptm.json" 
port="5000"

${CONSTRAINT_TOOLS}/constraint-tools dashboard \
  --model ${model} \
  --port ${port}

  
