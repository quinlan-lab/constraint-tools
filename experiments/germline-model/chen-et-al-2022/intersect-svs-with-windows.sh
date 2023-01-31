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

# Courtesy: Tom Nicholas 
# https://quinlangroup.slack.com/archives/D9Q0R5VNW/p1675111276956879?thread_ts=1674171985.879549&cid=D9Q0R5VNW
svs="/scratch/ucgd/lustre-work/quinlan/u0055382/SVAFotate/sources/SVAFotate_core_SV_popAFs.GRCh38.bed.gz"

windows="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/chen-mchale.bed"

svs_windows="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/svs-windows.bed"

get-svs () {
  zcat ${svs} \
    | python ${CONSTRAINT_TOOLS}/experiments/germline-model/chen-et-al-2022/filter_svs.py
}

get-svs-head () { 
  get-svs | head -1
}

get-svs-tail () { 
  get-svs | tail -n +2
}

get-windows-head () {
  head -1 ${windows}
}

get-windows-tail () {
  cat ${windows} | tail -n +2
}

header-line () {
  echo -e "$(get-svs-head)\t$(get-windows-head)"
}

intersect-svs-with-windows () {
  bedtools intersect \
    -a <(get-svs-tail) \
    -b <(get-windows-tail) \
    -wa \
    -wb 
}

(
  header-line 
  intersect-svs-with-windows
) > ${svs_windows}  

info "Wrote SVs with intersecting windows to:" ${svs_windows}  


