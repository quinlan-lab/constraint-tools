set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset 

source set-environment-variables.sh 

PATH="${CONSTRAINT_TOOLS}/bin:$PATH"

# https://stackoverflow.com/a/43476575/6674256
# need to export PYTHONPATH since it is not already in the environment: 
# `printenv | grep -w PYTHONPATH` returns zero output
export PYTHONPATH="${CONSTRAINT_TOOLS}/utilities"

get-merged-exons () {
  local merged_exons="${CONSTRAINT_TOOLS_DATA}/genes/grch38/exons.merged.bed"
  cat ${merged_exons} | uniq
}

get-exon-adjacent-windows () {
  chromosome_sizes="${CONSTRAINT_TOOLS_DATA}/reference/grch38/chromosome-sizes/hg38.chrom.sizes.sorted"
  bedtools flank -i <(get-merged-exons) -g ${chromosome_sizes} -b 1000 
}

get-trustworthy-exon-adjacent-windows () {
  trustworthy_noncoding_regions="${CONSTRAINT_TOOLS}/dist/trustworthy-noncoding-regions-germline-grch38.bed.gz"
  bedtools intersect -a <(get-exon-adjacent-windows) -b ${trustworthy_noncoding_regions} -f 1 -wa
}

output="${CONSTRAINT_TOOLS_DATA}/trustworthy-exon-adjacent-windows.bed"
get-trustworthy-exon-adjacent-windows > ${output}
info "Wrote trustworthy exon-adjacent windows to:" ${output}


