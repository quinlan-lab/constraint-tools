set -o errexit
set -o pipefail
set -o noclobber
set -o nounset

source set-environment-variables.sh 

model="${CONSTRAINT_TOOLS}/dist/model-germline-grch38.kmerSize-7.json" 

python ${CONSTRAINT_TOOLS}/utilities/read_model.py ${model}
