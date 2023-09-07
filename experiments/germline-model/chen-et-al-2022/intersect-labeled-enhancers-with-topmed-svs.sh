set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset 

source set-environment-variables.sh 

PATH="${CONSTRAINT_TOOLS}/experiments/germline-model/chen-et-al-2022:$PATH"

# https://stackoverflow.com/a/43476575/6674256
# need to export PYTHONPATH since it is not already in the environment: 
# `printenv | grep -w PYTHONPATH` returns zero output
export PYTHONPATH="${CONSTRAINT_TOOLS}/utilities"

LABELED_ENHANCERS="${CONSTRAINT_TOOLS_DATA}/khurana/labeled-enhancers.bed"
TOPMED_SVS="/scratch/ucgd/lustre-work/quinlan/u0055382/SVAFotate/supporting_data/TOPMed.GRCh38.bed.gz"

get-labeled-enhancers-head () {
  head -1 ${LABELED_ENHANCERS} 
}

get-labeled-enhancers-tail () {
  tail -n +2 ${LABELED_ENHANCERS} 
}

get-topmed-svs-head () {
  less ${TOPMED_SVS} \
    | head -1 
}

get-topmed-svs-tail () {
  less ${TOPMED_SVS} \
    | tail -n +2 \
    | awk '{print "chr"$0}'
}

intersect-enhancers-with-svs () {
  bedtools intersect \
    -a <(get-labeled-enhancers-tail) \
    -b <(get-topmed-svs-tail) \
    -wao  
}

create-header () {
  echo -e "$(get-labeled-enhancers-head)\t$(get-topmed-svs-head)\tenhancer-sv-overlap"
}

intersect-enhancers-with-svs-with-header () {
  enhancers_intersect_svs="${CONSTRAINT_TOOLS_DATA}/khurana/labeled-enhancers-intersect-topmed-svs.bed"
  (
    create-header
    intersect-enhancers-with-svs
  ) > ${enhancers_intersect_svs}  
  info "Wrote labeled enhancers with intersecting topmed svs to:" ${enhancers_intersect_svs}  
}

intersect-enhancers-with-svs-with-header
