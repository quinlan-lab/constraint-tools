# intersect Halldorsson windows with enhancers and exons 

# based on: 
# experiments/germline-model/chen-et-al-2022/add-overlapAmounts.sh

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

# Full path to the script
script_path=$(realpath "$0")

# Directory of the script
script_dir=$(dirname "$0")

halldorsson_windows="${CONSTRAINT_TOOLS_DATA}/depletion_rank_scores/41586_2022_4965_MOESM3_ESM"

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5467550/
# data courtesy: Tom Nicholas 
GeneHancer_enhancers="/scratch/ucgd/lustre-labs/quinlan/u0055382/genome_reference/genehancer/genehancer_GRCh38.bed.gz" 

merged_exons="${CONSTRAINT_TOOLS_DATA}/genes/grch38/exons.merged.bed"

get-halldorsson-windows-head () { 
  set +o errexit
  cat ${halldorsson_windows} | head -1 
  set -o errexit
}

get-halldorsson-windows-tail () { 
  cat ${halldorsson_windows} | tail -n +2
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

# get-halldorsson-windows-head 
# get-halldorsson-windows-tail | head 
# get-GeneHancer-enhancers | head
# get-merged-exons | head

augment-header-line () {
  get-halldorsson-windows-head | python ${script_dir}/augment_header_line.py 
}

# augment-header-line

add-overlapAmounts () {
  bedtools intersect \
      -a <(get-halldorsson-windows-tail) \
      -b \
        <(get-GeneHancer-enhancers) \
        <(get-merged-exons) \
      -names \
        enhancer \
        merged_exon \
      -wao 
}

# add-overlapAmounts | head 

halldorsson_windows_with_overlapAmounts="${CONSTRAINT_TOOLS_DATA}/depletion_rank_scores/41586_2022_4965_MOESM3_ESM.overlapAmounts.bed"
(
  augment-header-line 
  add-overlapAmounts
) > ${halldorsson_windows_with_overlapAmounts}
# | head -50 | column -t -s $'\t'   

info "Wrote overlap amounts to" ${halldorsson_windows_with_overlapAmounts}  

