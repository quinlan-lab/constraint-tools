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

OVERLAP_FRACTION=$1 

get-large-windows-head () {
  echo -e "chrom\tstart\tend"
}

get-large-windows-tail () {
  cat ${LARGE_WINDOWS} 
}

get-topmed-deletions-head () {
  less ${TOPMED_SVS} \
    | head -1 
}

# both het and homalt deletions: 
get-topmed-deletions-tail () {
  less ${TOPMED_SVS} \
    | tail -n +2 \
    | awk '{print "chr"$0}' \
    | awk '$5 == "DEL"' 
}

intersect-large-windows-with-deletions () {
  bedtools intersect \
    -a <(get-large-windows-tail) \
    -b <(get-topmed-deletions-tail) \
    -f ${OVERLAP_FRACTION} \
    -c
}

create-header () {
  echo -e "$(get-large-windows-head)\tnumber_of_overlapping_topmed_deletions"
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

write-large-windows-with-deletion-counts () {
  large_windows_intersect_deletions="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/filtered-large-windows-with-deletion-counts.bed"
  (
    create-header
    intersect-large-windows-with-deletions-and-filter
  ) > ${large_windows_intersect_deletions}
  info "Wrote (filtered) large windows with intersecting topmed deletions to:" ${large_windows_intersect_deletions}  
}  

write-large-windows-with-deletion-counts

