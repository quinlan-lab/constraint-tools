set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset

source set-environment-variables.sh 

PATH="${CONSTRAINT_TOOLS}/experiments/germline-model/chen-et-al-2022:$PATH" 

# https://stackoverflow.com/a/43476575/6674256
# need to export PYTHONPATH since it is not already in the environment: 
# `printenv | grep -w PYTHONPATH` returns zero output
export PYTHONPATH="${CONSTRAINT_TOOLS}/utilities"

exons="${CONSTRAINT_TOOLS_DATA}/genes/grch38/exons.sorted.bed.gz"
merged_exons="${CONSTRAINT_TOOLS_DATA}/genes/grch38/exons.merged.bed"

# zcat ${exons} | head 
# echo ""

merge-exons () {
  bedtools merge -i ${exons} 
}

merge-exons > ${merged_exons}
# | head 

info "Merged exons and wrote to:" ${merged_exons}
