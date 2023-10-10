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

LARGE_WINDOWS="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/large-windows.bed"
TOPMED_SVS="/scratch/ucgd/lustre-work/quinlan/u0055382/SVAFotate/supporting_data/TOPMed.GRCh38.bed.gz"

GET_TOPMED_DELETIONS_TAIL="${1}" 
DELETION_TYPE="${2}" 
LOWER_SIZE_LIMIT="${3}"
UPPER_SIZE_LIMIT="${4}"
ALLELE_FREQ_THRESHOLD="${5}"

get-large-windows-head () {
  echo -e "chrom_window\tstart_window\tend_window"
}

get-large-windows-tail () {
  cat ${LARGE_WINDOWS} 
}

get-topmed-deletions-head () {
  less ${TOPMED_SVS} \
    | head -1 
}

# both het and homalt deletions: 
get-all-topmed-deletions-tail () {
  less ${TOPMED_SVS} \
    | tail -n +2 \
    | awk '{print "chr"$0}' \
    | awk '$5 == "DEL"' 
}

get-rare-topmed-deletions-tail () { 
  get-all-topmed-deletions-tail \
    | awk \
      -v threshold=${ALLELE_FREQ_THRESHOLD} \
      '$8 < threshold'
}

get-common-topmed-deletions-tail () { 
  get-all-topmed-deletions-tail \
    | awk \
      -v threshold=${ALLELE_FREQ_THRESHOLD} \
      '$8 > threshold'
}

get-short-topmed-deletions-tail () { 
  get-all-topmed-deletions-tail \
    | awk \
      -v threshold=${LOWER_SIZE_LIMIT} \
      '$4 < threshold'
}

get-medium-topmed-deletions-tail () { 
  get-all-topmed-deletions-tail \
    | awk \
      -v lower=${LOWER_SIZE_LIMIT} \
      -v upper=${UPPER_SIZE_LIMIT} \
      '$4 >= lower && $4 < upper'
}

get-long-topmed-deletions-tail () { 
  get-all-topmed-deletions-tail \
    | awk \
      -v threshold=${UPPER_SIZE_LIMIT} \
      '$4 >= threshold'
}

intersect-large-windows-with-deletions () {
  bedtools intersect \
    -a <(get-large-windows-tail) \
    -b <(${GET_TOPMED_DELETIONS_TAIL}) \
    -c
}

create-header () {
  echo -e "$(get-large-windows-head)\tnumber_of_overlapping_topmed_deletions\tchrom_merged_deletion\tstart_merged_deletion\tend_merged_deletion\twindow_merged_deletion_overlap"
}

get-exclude-regions () {
  cat \
    "${CONSTRAINT_TOOLS_DATA}/chromosome-bands/grch38/centromeres.bed" \
    "${CONSTRAINT_TOOLS_DATA}/chromosome-bands/grch38/telomeres.bed" \
    <(zcat "${CONSTRAINT_TOOLS_DATA}/gaps/grch38/gaps.sorted.bed.gz") \
    | get-nonXY-chromosomes \
    | sort -k1,1 -k2,2n --version-sort
}

intersect-large-windows-with-deletions-and-filter () {
  bedtools subtract \
    -a <(intersect-large-windows-with-deletions) \
    -b <(get-exclude-regions) \
    -A 
}

get-merged-topmed-deletions-tail () {
  ${GET_TOPMED_DELETIONS_TAIL} \
    | bedtools merge -i - 
}

intersect-large-windows-with-merged-deletions () {
  bedtools intersect \
    -a <(intersect-large-windows-with-deletions-and-filter) \
    -b <(get-merged-topmed-deletions-tail) \
    -wao
}

write-large-windows-with-deletions () {
  local large_windows_intersect_deletions="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/filtered-large-windows-with-${DELETION_TYPE}-deletions.bed"
  (
    create-header
    intersect-large-windows-with-merged-deletions
  ) > ${large_windows_intersect_deletions}
  info "Wrote (filtered) large windows with intersecting ${DELETION_TYPE} topmed deletions to:" ${large_windows_intersect_deletions}  
}  

write-deletion-count () { 
  local number_deletions_filename="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/number-${DELETION_TYPE}-deletions.txt"
  ${GET_TOPMED_DELETIONS_TAIL} | wc -l > ${number_deletions_filename}
  info "Wrote number of deletions in this particular class to:" ${number_deletions_filename}
}

write-large-windows-with-deletions
write-deletion-count
