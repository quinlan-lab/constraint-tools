set -o errexit 
set -o nounset 
# set -o xtrace 

# https://stackoverflow.com/a/43476575/6674256
# need to export PYTHONPATH since it is not already in the environment: 
# `printenv | grep -w PYTHONPATH` returns zero output

CONSTRAINT_TOOLS=$1

export PYTHONPATH="${CONSTRAINT_TOOLS}/utilities:${CONSTRAINT_TOOLS}/predict-constraint"

# no need to export PATH since it is already in the environment: 
# `printenv | grep -w PATH` returns non-zero output 
PATH="${CONSTRAINT_TOOLS}/bin:$PATH" 

model="/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/model.ptm.json" 
region="chr1:100,000-100,600"
window_size="51" 
window_stride="25"

python ${CONSTRAINT_TOOLS}/predict-constraint/compute_mutation_counts.py \
    --region ${region} \
    --model ${model} \
    --window-size ${window_size} \
    --window-stride ${window_stride} \
  | jq .

