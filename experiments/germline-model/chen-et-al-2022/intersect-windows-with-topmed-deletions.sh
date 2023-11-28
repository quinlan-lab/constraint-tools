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

# bed file to intersect with SVs:
WINDOWS="${1}"

WINDOWS_WITH_DELETIONS="${2}" 
PROCESSED_DELETIONS="${3}"

TOPMED_SVS="/scratch/ucgd/lustre-work/quinlan/u0055382/SVAFotate/supporting_data/TOPMed.GRCh38.bed.gz"

# How to stratify deletions:
GET_TOPMED_DELETIONS_TAIL="${4}" 
DELETION_TYPE="${5}" 
LOWER_SIZE_LIMIT="${6}"
UPPER_SIZE_LIMIT="${7}"
ALLELE_FREQ_THRESHOLD="${8}"

# Filter out false deletions: 
SUSPICIOUS_DELETION_SIZE_THRESHOLD="1000000"

info "We assume that the first line of the following is a header line:" ${WINDOWS}

get-windows-head () {
  head -1 ${WINDOWS}
}

get-windows-tail () {
  tail -n +2 ${WINDOWS} 
}

get-topmed-deletions-head () {
  less ${TOPMED_SVS} \
    | head -1 
}

# both het and homalt deletions: 
# filter out suspiciously large deletions:
get-all-topmed-deletions-tail () {
  less ${TOPMED_SVS} \
    | tail -n +2 \
    | awk '{print "chr"$0}' \
    | awk '$5 == "DEL"' \
    | awk -v threshold=${SUSPICIOUS_DELETION_SIZE_THRESHOLD} '$4 < threshold'
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

get-windows-with-deletion-counts () {
  bedtools intersect \
    -a <(get-windows-tail) \
    -b <(${GET_TOPMED_DELETIONS_TAIL}) \
    -f 0.5 \
    -c
}

get-exclude-regions () {
  cat \
    "${CONSTRAINT_TOOLS_DATA}/chromosome-bands/grch38/centromeres.bed" \
    "${CONSTRAINT_TOOLS_DATA}/chromosome-bands/grch38/telomeres.bed" \
    <(zcat "${CONSTRAINT_TOOLS_DATA}/gaps/grch38/gaps.sorted.bed.gz") \
    | get-nonXY-chromosomes \
    | sort -k1,1 -k2,2n --version-sort
}

filter-windows-with-deletion-counts () {
  bedtools subtract \
    -a <(get-windows-with-deletion-counts) \
    -b <(get-exclude-regions) \
    -A 
}

get-merged-topmed-deletions-tail () {
  ${GET_TOPMED_DELETIONS_TAIL} \
    | bedtools merge -i - 
}

add-deletion-overlaps-to-windows () {
  bedtools intersect \
      -a <(filter-windows-with-deletion-counts) \
      -b <(get-merged-topmed-deletions-tail) \
      -wao \
    | cut -f1-7,11
}

create-header () {
  echo -e "$(get-windows-head)\tdeletion_count\tdeletion_overlap"
}

write-windows-with-deletions () {
  (
    create-header
    add-deletion-overlaps-to-windows
  ) > ${WINDOWS_WITH_DELETIONS}
  info "Wrote windows with deletion counts and overlaps to:" ${WINDOWS_WITH_DELETIONS}  
}  

write-deletion-count () { 
  local number_deletions_filename="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/number-${DELETION_TYPE}-deletions.txt"
  ${GET_TOPMED_DELETIONS_TAIL} | wc -l > ${number_deletions_filename}
  info "Wrote number of deletions in this particular stratum to:" ${number_deletions_filename}
}

write-deletions () { 
  ${GET_TOPMED_DELETIONS_TAIL} | cut -f1-3 > ${PROCESSED_DELETIONS}
  info "Wrote deletions in this particular stratum to:" ${PROCESSED_DELETIONS}
}

write-windows-with-deletions
write-deletion-count
write-deletions
