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

chen_mchale_windows="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/chen-mchale.multiple-kmers.bed"

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5467550/
# data courtesy: Tom Nicholas 
GeneHancer_enhancers="/scratch/ucgd/lustre-work/quinlan/u0055382/genome_reference/genehancer/genehancer_GRCh38.bed.gz" 

merged_exons="${CONSTRAINT_TOOLS_DATA}/genes/grch38/exons.merged.bed"
gaps="${CONSTRAINT_TOOLS_DATA}/gaps/grch38/gaps.sorted.bed.gz"
encode_exclude_regions="${CONSTRAINT_TOOLS_DATA}/encode-exclude-regions/grch38/encode-exclude-regions.sorted.bed.gz"
covered_regions="${CONSTRAINT_TOOLS_DATA}/gnomad/v3/gnomad_v3_coverage.filtered.sorted.bed.gz"

chromosome_sizes="${CONSTRAINT_TOOLS_DATA}/reference/grch38/chromosome-sizes/hg38.chrom.sizes.sorted"

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

get-merged-gaps () {
  bedtools merge -i ${gaps} 
}

get-merged-encode-exclude-regions () {
  bedtools merge -i ${encode_exclude_regions} 
}

get-uncovered-regions () {
  bedtools merge -i ${covered_regions} \
    | bedtools complement -i - -g ${chromosome_sizes} \
    | awk '$3 - $2 > 50'
}

# get-chen-windows | head 
# get-GeneHancer-enhancers | head
# get-merged-exons | head
# get-merged-gaps | head
# get-merged-encode-exclude-regions | head
# get-uncovered-regions | wc -l 
# exit 1 

header-line () {
  python experiments/germline-model/chen-et-al-2022/augment_header_line.py  
}

add-overlapAmounts () {
  bedtools intersect \
      -a <(get-chen-windows) \
      -b \
        <(get-GeneHancer-enhancers) \
        <(get-merged-exons) \
        <(get-merged-gaps) \
        <(get-merged-encode-exclude-regions) \
        <(get-uncovered-regions) \
      -names \
        enhancer \
        merged_exon \
        merged_gap \
        merged_encode_exclude_region \
        uncovered_region \
      -wao 
}

chen_mchale_windows_with_overlapAmounts="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/chen-mchale.multiple-kmers.overlapAmounts.bed"
(
  header-line 
  add-overlapAmounts
) > ${chen_mchale_windows_with_overlapAmounts} 
# | head -50 | column -t -s $'\t' 

info "Wrote overlap amounts to" ${chen_mchale_windows_with_overlapAmounts}  

