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

get-svs () {
  local source=$1
  local number_individuals_who_were_genotyped_min=$2
  # Courtesy: Tom Nicholas 
  # https://quinlangroup.slack.com/archives/D9Q0R5VNW/p1675111276956879?thread_ts=1674171985.879549&cid=D9Q0R5VNW
  local svs="/scratch/ucgd/lustre-work/quinlan/u0055382/SVAFotate/sources/SVAFotate_core_SV_popAFs.GRCh38.bed.gz"
  zcat ${svs} \
    | python ${CONSTRAINT_TOOLS}/experiments/germline-model/chen-et-al-2022/filter_svs.py ${source} ${number_individuals_who_were_genotyped_min}
}

get-svs-head () { 
  local source=$1
  local number_individuals_who_were_genotyped_min=$2
  get-svs ${source} ${number_individuals_who_were_genotyped_min} | head -1
}

get-svs-tail () { 
  local source=$1
  local number_individuals_who_were_genotyped_min=$2
  get-svs ${source} ${number_individuals_who_were_genotyped_min} | tail -n +2
}

get-windows () {                                                                                            
  local windows="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/chen-mchale.kmerSizes.trainSets.enhancer-exon.bed"
  cat ${windows}
}

get-windows-head () {
  get-windows | head -1
}

get-windows-tail () {
  get-windows | tail -n +2
}

header-line () {
  local source=$1
  local number_individuals_who_were_genotyped_min=$2
  echo -e "$(get-svs-head ${source} ${number_individuals_who_were_genotyped_min})\t$(get-windows-head)"
}

get-merged-exons () {
  local merged_exons="${CONSTRAINT_TOOLS_DATA}/genes/grch38/exons.merged.bed"
  cat ${merged_exons} | uniq
}

get-noncoding-svs-tail () { 
  local source=$1
  local number_individuals_who_were_genotyped_min=$2
  bedtools intersect \
    -a <(get-svs-tail ${source} ${number_individuals_who_were_genotyped_min}) \
    -b <(get-merged-exons) \
    -v \
    -wa 
}

intersect-noncoding-svs-with-windows () {
  local source=$1
  local number_individuals_who_were_genotyped_min=$2
  bedtools intersect \
    -a <(get-noncoding-svs-tail ${source} ${number_individuals_who_were_genotyped_min}) \
    -b <(get-windows-tail) \
    -wa \
    -wb 
}

intersect-noncoding-svs-with-windows-with-header () {
  local source=$1
  local number_individuals_who_were_genotyped_min=$2
  noncoding_svs_windows="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/${source}-noncoding-svs-chen-mchale.kmerSizes.trainSets.enhancer-exon.bed"
  (
    header-line ${source} ${number_individuals_who_were_genotyped_min}
    intersect-noncoding-svs-with-windows ${source} ${number_individuals_who_were_genotyped_min}
  ) > ${noncoding_svs_windows}  
  info "Wrote ${source} noncoding SVs with intersecting windows to:" ${noncoding_svs_windows}  
}

intersect-noncoding-svs-with-windows-with-header "gnomAD" "10000"
intersect-noncoding-svs-with-windows-with-header "1000G" "2500"
intersect-noncoding-svs-with-windows-with-header "CCDG" "n/a"



