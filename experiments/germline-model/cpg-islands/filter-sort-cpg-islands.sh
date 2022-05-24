set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset 

source download-process-data/set-environment-variables.sh 

# no need to export PATH since it is already in the environment: 
# `printenv | grep -w PATH` returns non-zero output 
PATH="${CONSTRAINT_TOOLS}/bin:$PATH" 

info "jq path is:" "$(which jq)"
info "bedtools path is:" $(which bedtools)

model="${CONSTRAINT_TOOLS}/dist/model-germline-grch38.json"
neutral_regions=$(jq --raw-output .neutralRegions ${model})
exons="${CONSTRAINT_TOOLS_DATA}/genes/grch38/canonical-exons.sorted.bed.gz"
CpGs="${CONSTRAINT_TOOLS_DATA}/cpg-islands/grch38/cpg-islands.sorted.bed.gz"

CpG_island_length_column="5"

percentage_of_island_that_is_CpG_column="8"
percentage_of_island_that_is_C_or_G_column="9"
observed_CpG_number_relative_to_expectation_column="10"

filter () {
  zcat ${CpGs} |
    awk -F $'\t' -v CpG_island_length=${CpG_island_length_column} '$CpG_island_length > 1000' |
    bedtools intersect -a - -b ${neutral_regions} -wa -f 0.5 -u
}

filter_and_sort () {
  local sort_column=$1
  filter | sort --field-separator=$'\t' -k${sort_column},${sort_column}n
}

filter_and_sort ${percentage_of_island_that_is_CpG_column}

