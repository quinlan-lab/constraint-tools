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

WINDOWS="${1}"
DELETIONS_CLASS="${2}"
WINDOWS_WITH_DELETIONS="${3}" 

PUBLIC_REPO_DIR="${4}"
PROCESSED_DELETIONS="${5}"

# How to stratify deletions:
DELETION_TYPE="${6}" 
GET_DELETIONS_TAIL="get-${DELETION_TYPE}-deletions-tail"

LOWER_SIZE_LIMIT="${7}"
UPPER_SIZE_LIMIT="${8}"
ALLELE_FREQ_THRESHOLD="${9}"

# Filter out false deletions: 
SUSPICIOUS_DELETION_SIZE_THRESHOLD="1000000"

info "We assume that the first line of the following is a header line:" ${WINDOWS}

source "${CONSTRAINT_TOOLS}/experiments/germline-model/chen-et-al-2022/get-${DELETIONS_CLASS}-deletions.sh"

get-windows-head () {
  head -1 ${WINDOWS}
}

get-windows-tail () {
  tail -n +2 ${WINDOWS} 
}

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

get-windows-with-deletion-counts () {
  bedtools intersect \
    -a <(get-windows-tail) \
    -b <(${GET_DELETIONS_TAIL}) \
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

get-merged-deletions-tail () {
  ${GET_DELETIONS_TAIL} \
    | bedtools merge -i - 
}

add-deletion-overlaps-to-windows () {
  bedtools intersect \
      -a <(filter-windows-with-deletion-counts) \
      -b <(get-merged-deletions-tail) \
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

push-deletions-to-public-repo () { 
  local output="${PUBLIC_REPO_DIR}/${PROCESSED_DELETIONS}"
  echo "track name=${DELETION_TYPE}-${DELETIONS_CLASS}-deletions description=${DELETION_TYPE}-${DELETIONS_CLASS}-deletions color=255,0,0," > ${output}
  ${GET_DELETIONS_TAIL} | cut -f1-3 >> ${output}
  info "Wrote ${DELETION_TYPE} ${DELETIONS_CLASS} deletions in UCSC-genome-browser format to:" ${output}

  cd ${PUBLIC_REPO_DIR}
  git add ${output}
  ( git commit -m "Add ${DELETION_TYPE} ${DELETIONS_CLASS} deletions" ) || true
  git push
  info "Pushed ${output} to public repo"
}

write-windows-with-deletions
push-deletions-to-public-repo
