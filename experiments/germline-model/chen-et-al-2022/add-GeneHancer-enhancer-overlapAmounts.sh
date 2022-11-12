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

get-chen-windows () { 
  cat ${chen_mchale_windows} | tail -n +2
}

get-GeneHancer-enhancers () {
  zcat ${GeneHancer_enhancers} \
    | tail -n +2 \
    | cut -f1-3 \
    | awk -v OFS="\t" '{ print "chr"$1, $2, $3}' \
    | awk -v min_enhancer_size=${min_enhancer_size} '$3-$2 > min_enhancer_size' \
    | awk -v max_enhancer_size=${max_enhancer_size} '$3-$2 < max_enhancer_size' \
    | sort --version-sort -k1,1 -k2,2n \
    | uniq 
}

# get-chen-windows | head 
# min_enhancer_size=0
# max_enhancer_size=750
# get-GeneHancer-enhancers | head
# exit 1 

header-line () {
  echo -e "chromosome\tstart\tend\tposition\tN_bar\tN_observed\tK_bar\tK_observed\tM\tchen_zscore\tenhancer_chromosome\tenhancer_start\tenhancer_end\tchen_enhancer_overlap_bps"
}

add-GeneHancer-enhancer-overlapAmounts () {
  bedtools intersect \
      -a <(get-chen-windows) \
      -b <(get-GeneHancer-enhancers) \
      -wao 
}

for enhancer_size_bounds in "0,750" "750,100000" "0,100000"; do
  chen_mchale_windows_with_GeneHancer_enhancer_overlapAmounts="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/chen-mchale.GeneHancer-enhancers.${enhancer_size_bounds}.overlapAmounts.bed"
  IFS="," read min_enhancer_size max_enhancer_size <<< ${enhancer_size_bounds}
  (
    header-line 
    add-GeneHancer-enhancer-overlapAmounts
  ) > ${chen_mchale_windows_with_GeneHancer_enhancer_overlapAmounts}   
  # | head | column -t -s $'\t' 
  info "Wrote GeneHancer-enhancer overlap amounts to" ${chen_mchale_windows_with_GeneHancer_enhancer_overlapAmounts}  
done

