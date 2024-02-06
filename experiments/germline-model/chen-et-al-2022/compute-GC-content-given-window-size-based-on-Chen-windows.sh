set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset

source /scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools/set-environment-variables.sh 

PATH="${CONSTRAINT_TOOLS}/experiments/germline-model/chen-et-al-2022:$PATH"

# https://stackoverflow.com/a/43476575/6674256
# need to export PYTHONPATH since it is not already in the environment: 
# `printenv | grep -w PYTHONPATH` returns zero output
export PYTHONPATH="${CONSTRAINT_TOOLS}/utilities"

CHEN_DATA_DIRECTORY="${CONSTRAINT_TOOLS_DATA}/chen-et-al-2023-published-version/41586_2023_6045_MOESM4_ESM"
CHEN_FILE_STEM="Supplementary_Data_2" 

WINDOW_SIZE="${1}" 

source "${CONSTRAINT_TOOLS}/experiments/germline-model/chen-et-al-2022/get-windows-from-chen.sh"

compute-GC-content-head () {
  get-windows-head \
    | awk '{ print $0"\twindow_GC_content" }'
}

compute-GC-content-tail () {
  local genome="${CONSTRAINT_TOOLS_DATA}/reference/grch38/hg38.analysisSet.fa"
  bedtools nuc \
      -fi ${genome} \
      -bed <(get-windows-tail) \
    | cut -f1-7,9 \
    | tail -n +2
}

write-to-disk () {
  local filename="${CHEN_DATA_DIRECTORY}/${CHEN_FILE_STEM}.GC_content_${WINDOW_SIZE}.bed"
  # set +o errexit
  (
    compute-GC-content-head
    compute-GC-content-tail 
  ) > ${filename}
  # set -o errexit
  info "Wrote GC content:" ${filename}
}

write-to-disk 
