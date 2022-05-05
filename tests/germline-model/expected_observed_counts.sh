set -o errexit 
set -o nounset 
# set -o xtrace 

# https://stackoverflow.com/a/43476575/6674256
# need to export PYTHONPATH since it is not already in the environment: 
# `printenv | grep -w PYTHONPATH` returns zero output

CONSTRAINT_TOOLS="/scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools"

export PYTHONPATH="${CONSTRAINT_TOOLS}/utilities:${CONSTRAINT_TOOLS}/predict-constraint"

# no need to export PATH since it is already in the environment: 
# `printenv | grep -w PATH` returns non-zero output 
PATH="${CONSTRAINT_TOOLS}/bin:$PATH" 

model="/scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools/dist/model-germline-grch38.json" 
region="chr1:1,000,000-1,000,200"
window_stride="25"

python ${CONSTRAINT_TOOLS}/predict-constraint/germline-model/expected_observed_counts.py \
    --region ${region} \
    --model ${model} \
    --window-stride ${window_stride} 
