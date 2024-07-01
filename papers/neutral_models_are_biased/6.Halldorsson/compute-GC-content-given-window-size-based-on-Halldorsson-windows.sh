set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset

source /scratch/ucgd/lustre-labs/quinlan/u6018199/constraint-tools/set-environment-variables.sh 

# https://stackoverflow.com/a/43476575/6674256
# need to export PYTHONPATH since it is not already in the environment: 
# `printenv | grep -w PYTHONPATH` returns zero output
export PYTHONPATH="${CONSTRAINT_TOOLS}/utilities"

HALLDORSSON_DATA_DIRECTORY="${CONSTRAINT_TOOLS_DATA}/depletion_rank_scores"
HALLDORSSON_FILE_STEM="41586_2022_4965_MOESM3_ESM" 

WINDOW_SIZE="${1}" 

source "${CONSTRAINT_TOOLS}/papers/neutral_models_are_biased/6.Halldorsson/get-windows-from-halldorsson.sh"

compute-GC-content-head () {
  get-windows-head \
    | awk '{ print $0"\twindow_GC_content" }'
}

# "tail" is required because "bedtools nuc" generates a header 
compute-GC-content-tail () {
  local genome="${CONSTRAINT_TOOLS_DATA}/reference/grch38/hg38.analysisSet.fa"
  bedtools nuc \
      -fi ${genome} \
      -bed <(get-windows-tail) \
    | cut -f1-7,9 \
    | tail -n +2
}

write-to-disk () {
  local filename="${HALLDORSSON_DATA_DIRECTORY}/${HALLDORSSON_FILE_STEM}.GC_content_${WINDOW_SIZE}.bed"
  info "Writing GC content to:" ${filename}
  # set +o errexit
  (
    compute-GC-content-head
    compute-GC-content-tail 
  ) > ${filename}
  # set -o errexit
  info "Wrote GC content:" ${filename}
}

write-to-disk 
