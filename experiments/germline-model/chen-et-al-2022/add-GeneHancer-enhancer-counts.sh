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

chen_mchale_windows_with_GeneHancer_enhancer_counts="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/chen-mchale.GeneHancer-enhancers.bed"

header-line () {
  echo -e "chromosome\tchen_start\tchen_end\tchen_zscore\tmchale_start\tmchale_end\tmchale_position\tmchale_N_bar\tmchale_N_observed\tmchale_K_bar\tmchale_K_observed\tmchale_M\toverlap_bps\tchen_GeneHancer_enhancer_count"
}

get-chen-windows () { 
  cat ${chen_mchale_windows} | tail -n +2
}

get-GeneHancer-enhancers () {
  zcat ${GeneHancer_enhancers} \
    | tail -n +2 \
    | cut -f1-3 \
    | awk -v OFS="\t" '{ print "chr"$1, $2, $3 }' \
    | sort --version-sort -k1,1 -k2,2n \
    | uniq 
}

# get-chen-windows | head 
# get-GeneHancer-enhancers | head 

add-GeneHancer-enhancer-counts () {
  bedtools intersect \
      -a <(get-chen-windows) \
      -b <(get-GeneHancer-enhancers) \
      -c
      # -wo 
}

(
  header-line 
  add-GeneHancer-enhancer-counts
) > ${chen_mchale_windows_with_GeneHancer_enhancer_counts}  
# | head | column -t -s $'\t'

info "Wrote GeneHancer enhancer counts to" ${chen_mchale_windows_with_GeneHancer_enhancer_counts}  


