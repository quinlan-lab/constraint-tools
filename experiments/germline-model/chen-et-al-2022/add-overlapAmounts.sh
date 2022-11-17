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

chen_mchale_windows="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/chen-mchale.bed"

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5467550/
# data courtesy: Tom Nicholas 
GeneHancer_enhancers="/scratch/ucgd/lustre-work/quinlan/u0055382/genome_reference/genehancer/genehancer_GRCh38.bed.gz" 

merged_exons="${CONSTRAINT_TOOLS_DATA}/genes/grch38/exons.merged.bed"

get-chen-windows () { 
  cat ${chen_mchale_windows} | tail -n +2
}

get-GeneHancer-enhancers () {
  zcat ${GeneHancer_enhancers} \
    | tail -n +2 \
    | cut -f1-3 \
    | awk -v OFS="\t" '{ print "chr"$1, $2, $3}' \
    | sort --version-sort -k1,1 -k2,2n \
    | uniq 
}

get-merged-exons () {
  cat ${merged_exons} | uniq 
}

# get-chen-windows | head 
# get-GeneHancer-enhancers | head
# get-merged-exons | head
# exit 1 

header-line () {
  echo -e "chromosome\tstart\tend\tposition\tN_bar\tN_observed\tK_bar\tK_observed\tM\tchen_zscore\tfeature\tfeature_chromosome\tfeature_start\tfeature_end\twindow_feature_overlap_bps"
}

add-overlapAmounts () {
  bedtools intersect \
      -a <(get-chen-windows) \
      -b <(get-GeneHancer-enhancers) <(get-merged-exons) \
      -names enhancer merged_exon \
      -wao 
}

chen_mchale_windows_with_overlapAmounts="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/chen-mchale.overlapAmounts.bed"
(
  header-line 
  add-overlapAmounts
) > ${chen_mchale_windows_with_overlapAmounts}
# | head -50 | column -t -s $'\t' 

info "Wrote overlap amounts to" ${chen_mchale_windows_with_overlapAmounts}  

