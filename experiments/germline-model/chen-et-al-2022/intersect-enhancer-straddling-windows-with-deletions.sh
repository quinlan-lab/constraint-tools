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

DELETIONS_CLASS="${1}"
WINDOWS_WITH_DELETIONS="${2}" 

# How to stratify deletions:
DELETION_TYPE="${3}" 
GET_DELETIONS_TAIL="get-${DELETION_TYPE}-deletions-tail"

LOWER_SIZE_LIMIT="${4}"
UPPER_SIZE_LIMIT="${5}"
ALLELE_FREQ_THRESHOLD="${6}"

# Filter out false deletions: 
SUSPICIOUS_DELETION_SIZE_THRESHOLD="${7}"

WINDOW_SIZE="${8}" 

ENHANCERS_CLASS="${9}"

CHROMOSOME_SIZES="${CONSTRAINT_TOOLS_DATA}/reference/grch38/chromosome-sizes/hg38.chrom.sizes.sorted"

source "${CONSTRAINT_TOOLS}/experiments/germline-model/chen-et-al-2022/get-${DELETIONS_CLASS}-deletions.sh"

source "${CONSTRAINT_TOOLS}/experiments/germline-model/chen-et-al-2022/get-windows-from-${ENHANCERS_CLASS}-enhancers.sh"

get-all-deletions-tail () { 
  get-deletions-tail
}

get-rare-deletions-tail () { 
  get-deletions-tail \
    | awk \
      -v threshold=${ALLELE_FREQ_THRESHOLD} \
      '$8 < threshold'
}

get-common-deletions-tail () { 
  get-deletions-tail \
    | awk \
      -v threshold=${ALLELE_FREQ_THRESHOLD} \
      '$8 > threshold'
}

get-short-deletions-tail () { 
  get-deletions-tail \
    | awk \
      -v threshold=${LOWER_SIZE_LIMIT} \
      '$4 < threshold'
}

get-medium-deletions-tail () { 
  get-deletions-tail \
    | awk \
      -v lower=${LOWER_SIZE_LIMIT} \
      -v upper=${UPPER_SIZE_LIMIT} \
      '$4 >= lower && $4 < upper'
}

get-long-deletions-tail () { 
  get-deletions-tail \
    | awk \
      -v threshold=${UPPER_SIZE_LIMIT} \
      '$4 >= threshold'
}

get-windows-with-deletion-overlaps-head () {
  echo -e "$(get-windows-head)\t$(get-deletions-head)"
}  

get-windows-with-deletion-overlaps-tail () {
  bedtools intersect \
    -a <(get-windows-tail) \
    -b <(${GET_DELETIONS_TAIL}) \
    -wa -wb
}

get-exclude-regions () {
  cat \
    "${CONSTRAINT_TOOLS_DATA}/chromosome-bands/grch38/centromeres.bed" \
    "${CONSTRAINT_TOOLS_DATA}/chromosome-bands/grch38/telomeres.bed" \
    <(zcat "${CONSTRAINT_TOOLS_DATA}/gaps/grch38/gaps.sorted.bed.gz") \
    | get-nonXY-chromosomes \
    | sort -k1,1 -k2,2n --version-sort
}

filter-windows-with-deletion-overlaps-tail () {
  bedtools subtract \
    -a <(get-windows-with-deletion-overlaps-tail) \
    -b <(get-exclude-regions) \
    -A 
}

write-windows-with-deletion-overlaps () {
  (
    get-windows-with-deletion-overlaps-head
    filter-windows-with-deletion-overlaps-tail
  ) > ${WINDOWS_WITH_DELETIONS}
  info "Wrote windows with deletion overlaps to:" ${WINDOWS_WITH_DELETIONS}  
}  

write-windows-with-deletion-overlaps
