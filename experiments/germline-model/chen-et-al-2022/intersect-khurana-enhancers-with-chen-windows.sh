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

get-khurana-enhancers () {
  local khurana_enhancers="${CONSTRAINT_TOOLS_DATA}/khurana/enhancers-and-khurana-scores.hg38.sorted.bed"
  cat ${khurana_enhancers} 
}

get-chen-windows () {                   
  cat ${CONSTRAINT_TOOLS_DATA}/chen-et-al-2022-second-preprint/Supplementary_Datasets/Supplementary_Data_2.bed 
}

intersect-enhancers-with-windows () {
  khurana_enhancers_intersect_chen_windows="${CONSTRAINT_TOOLS_DATA}/khurana/khurana-enhancers-intersect-chen-windows.bed"
  bedtools intersect \
      -a <(get-khurana-enhancers) \
      -b <(get-chen-windows) \
      -wa \
      -wb \
    > ${khurana_enhancers_intersect_chen_windows}  
  info "Wrote khurana enhancers with intersecting chen windows to:" ${khurana_enhancers_intersect_chen_windows}  
}

intersect-enhancers-with-windows



