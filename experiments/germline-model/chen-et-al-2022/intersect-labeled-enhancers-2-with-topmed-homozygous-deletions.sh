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

LABELED_ENHANCERS="${CONSTRAINT_TOOLS_DATA}/khurana/labeled-enhancers.2.bed"
TOPMED_SVS="/scratch/ucgd/lustre-work/quinlan/u0055382/SVAFotate/supporting_data/TOPMed.GRCh38.bed.gz"

get-labeled-enhancers-head () {
  head -1 ${LABELED_ENHANCERS} 
}

get-labeled-enhancers-tail () {
  tail -n +2 ${LABELED_ENHANCERS} 
}

get-topmed-homozygous-deletions-head () {
  less ${TOPMED_SVS} \
    | head -1 
}

get-topmed-homozygous-deletions-tail () {
  less ${TOPMED_SVS} \
    | tail -n +2 \
    | awk '{print "chr"$0}' \
    | awk '$5 == "DEL" && $11 > 0' 
}

intersect-enhancers-with-deletions () {
  bedtools intersect \
    -a <(get-labeled-enhancers-tail) \
    -b <(get-topmed-homozygous-deletions-tail) \
    -f 1.0 \
    -wao  
}

create-header () {
  echo -e "$(get-labeled-enhancers-head)\t$(get-topmed-homozygous-deletions-head)\tenhancer-deletion-overlap"
}

intersect-enhancers-with-deletions-with-header () {
  enhancers_intersect_deletions="${CONSTRAINT_TOOLS_DATA}/khurana/labeled-enhancers-2-intersect-topmed-homozygous-deletions.bed"
  (
    create-header
    intersect-enhancers-with-deletions
  ) > ${enhancers_intersect_deletions}  
  info "Wrote labeled enhancers with intersecting topmed homozygous deletions to:" ${enhancers_intersect_deletions}  
}

intersect-enhancers-with-deletions-with-header
