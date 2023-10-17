set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset 

source /scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools/set-environment-variables.sh 

PATH="${CONSTRAINT_TOOLS}/experiments/germline-model/chen-et-al-2022:$PATH"

LABELED_REGIONS_VERSION=$1
DELETION_SIZE_THRESHOLD=$2

# https://stackoverflow.com/a/43476575/6674256
# need to export PYTHONPATH since it is not already in the environment: 
# `printenv | grep -w PYTHONPATH` returns zero output
export PYTHONPATH="${CONSTRAINT_TOOLS}/utilities"

LABELED_REGIONS="${CONSTRAINT_TOOLS_DATA}/khurana/labeled-regions.${LABELED_REGIONS_VERSION}.bed"
TOPMED_SVS="/scratch/ucgd/lustre-work/quinlan/u0055382/SVAFotate/supporting_data/TOPMed.GRCh38.bed.gz"

get-labeled-regions-head () {
  head -1 ${LABELED_REGIONS} 
}

get-labeled-regions-tail () {
  tail -n +2 ${LABELED_REGIONS} 
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
    | awk '$5 == "DEL"' \
    | awk -v threshold=${DELETION_SIZE_THRESHOLD} '$4 < threshold'
}

intersect-regions-with-deletions () {
  bedtools intersect \
    -a <(get-labeled-regions-tail) \
    -b <(get-topmed-deletions-tail) \
    -f 1.0 \
    -wao  
}

create-header () {
  echo -e "$(get-labeled-regions-head)\t$(get-topmed-deletions-head)\tregion-deletion-overlap"
}

intersect-regions-with-deletions-with-header () {
  regions_intersect_deletions="${CONSTRAINT_TOOLS_DATA}/khurana/labeled-regions-${LABELED_REGIONS_VERSION}-intersect-sizeFiltered-topmed-deletions.bed"
  (
    create-header
    intersect-regions-with-deletions
  ) > ${regions_intersect_deletions}  
  info "Wrote labeled regions with intersecting size-filtered topmed deletions to:" ${regions_intersect_deletions}  
}

intersect-regions-with-deletions-with-header
