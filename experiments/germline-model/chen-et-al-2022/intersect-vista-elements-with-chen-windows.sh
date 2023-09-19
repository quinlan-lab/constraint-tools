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

get-vista-elements () {
  local label=$1
  # download-process-data/vista-elements/
  local vista_elements="${CONSTRAINT_TOOLS_DATA}/vista-elements/vista-elements.${label}.hg38.hg19.tsv"
  cut -f1-3 ${vista_elements} 
}

get-chen-windows () {                   
  cat ${CONSTRAINT_TOOLS_DATA}/chen-et-al-2022-second-preprint/Supplementary_Datasets/Supplementary_Data_2.bed 
}

intersect-vista-elements-with-chen-windows () {
  local label=$1
  vista_elements_intersect_chen_windows="${CONSTRAINT_TOOLS_DATA}/vista-elements/vista-elements-${label}-intersect-chen-windows.bed"
  bedtools intersect \
      -a <(get-vista-elements ${label}) \
      -b <(get-chen-windows) \
      -wa \
      -wb \
    > ${vista_elements_intersect_chen_windows}  
  info "Wrote vista-${label} elements with intersecting chen windows to:" ${vista_elements_intersect_chen_windows}  
}

intersect-vista-elements-with-chen-windows "positive" 
intersect-vista-elements-with-chen-windows "negative"




