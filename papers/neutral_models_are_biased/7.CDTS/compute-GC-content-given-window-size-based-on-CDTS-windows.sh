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

CDTS_DATA_DIRECTORY="${CONSTRAINT_TOOLS_DATA}/CDTS"
CDTS_FILE_STEM="CDTS.gnomAD.hg38.noncoding" 

WINDOW_SIZE="${1}" 

source "${CONSTRAINT_TOOLS}/papers/neutral_models_are_biased/7.CDTS/get-windows-from-CDTS.sh"

compute-GC-content-head () {
  get-GC-windows-head \
    | awk '{ print $0"\tGC_window__GC_content" }'
}

# "tail" is required because "bedtools nuc" generates a header 
compute-GC-content-tail () {
  local genome="${CONSTRAINT_TOOLS_DATA}/reference/grch38/hg38.analysisSet.fa"
  bedtools nuc \
      -fi ${genome} \
      -bed <(get-GC-windows-tail) \
    | cut -f1-6,8 \
    | tail -n +2
}

write-to-disk () {
  local filename="${CDTS_DATA_DIRECTORY}/${CDTS_FILE_STEM}.GC_content_${WINDOW_SIZE}.bed"
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
