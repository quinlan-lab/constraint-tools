set -o errexit
set -o pipefail
set -o noclobber
set -o nounset

tmpdir="/scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools/tests/germline-model/tmpdir"
model="/scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools/tests/germline-model/model-germline-grch38.json"
progress_bars='disk'
kmer_size="3"

mkdir --parents ${tmpdir}

source set-environment-variables.sh 

PATH="${CONSTRAINT_TOOLS}/train/germline-model:$PATH" 

# https://stackoverflow.com/a/43476575/6674256
# need to export PYTHONPATH since it is not already in the environment: 
# `printenv | grep -w PYTHONPATH` returns zero output
export PYTHONPATH="${CONSTRAINT_TOOLS}/utilities"

estimate-probabilities \
  --genome "XXX" \
  --mutations "XXX" \
  --number-chromosomes-min "-1" \
  --kmer-size ${kmer_size} \
  --tmpdir ${tmpdir} \
  --window-size "-1" \
  --model ${model} \
  --neutral-regions "XXX" \
  --progress-bars "${progress_bars}"

