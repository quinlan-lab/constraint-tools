# intersect CDTS windows with enhancers and exons 

# based on: 
# experiments/germline-model/chen-et-al-2022/add-overlapAmounts.sh

set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset

source set-environment-variables.sh 

# https://stackoverflow.com/a/43476575/6674256
# need to export PYTHONPATH since it is not already in the environment: 
# `printenv | grep -w PYTHONPATH` returns zero output
export PYTHONPATH="${CONSTRAINT_TOOLS}/utilities"

# Directory of the script
script_dir=$(dirname "$0")

CDTS_windows="${CONSTRAINT_TOOLS_DATA}/CDTS/CDTS.gnomAD.hg38.header.bed"

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5467550/
# data courtesy: Tom Nicholas 
GeneHancer_enhancers="/scratch/ucgd/lustre-labs/quinlan/u0055382/genome_reference/genehancer/genehancer_GRCh38.bed.gz" 

merged_exons="${CONSTRAINT_TOOLS_DATA}/genes/grch38/exons.merged.bed"

chromosome_sizes="${CONSTRAINT_TOOLS_DATA}/reference/grch38/chromosome-sizes/hg38.chrom.sizes.sorted"

get-CDTS-windows-head () { 
  set +o errexit
  cat ${CDTS_windows} | head -1 
  set -o errexit
}

get-CDTS-windows-tail () { 
  # "the calculated CDTS across the 550-bp window was attributed to the middle 10-bp bin"
  # omit alternate chromosomes, and X and Y chromosomes 
  cat ${CDTS_windows} \
    | tail -n +2 \
    | awk '{if ($3 - $2 == 10) print $0}' \
    | awk  \
      --assign OFS=$'\t' \
      '{ 
        midpoint = 0.5*($2+$3)
        start = midpoint-1 
        end = midpoint
        printf "%s\t%d\t%d", $1, start, end
        for(i=4; i<=NF; i++) {
          printf "\t%s", $i
        }
        printf "\n"
      }' \
    | get-nonXY-chromosomes \
    | bedtools slop -i - -g ${chromosome_sizes} -b 275
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

augment-header-line () {
  get-CDTS-windows-head | python ${script_dir}/augment_header_line.py 
}

add-overlapAmounts () {
  bedtools intersect \
      -a <(get-CDTS-windows-tail) \
      -b \
        <(get-GeneHancer-enhancers) \
        <(get-merged-exons) \
      -names \
        enhancer \
        merged_exon \
      -wao 
}

CDTS_windows_with_overlapAmounts="${CONSTRAINT_TOOLS_DATA}/CDTS/CDTS.gnomAD.hg38.overlapAmounts.bed"
(
  augment-header-line 
  add-overlapAmounts
) > ${CDTS_windows_with_overlapAmounts}
# | head -50 | column -t -s $'\t'   

info "Wrote overlap amounts to" ${CDTS_windows_with_overlapAmounts}  

