set -o errexit
set -o pipefail
# set -o xtrace
set -o nounset

source download-data/set-environment-variables.sh 

exons="${CONSTRAINT_TOOLS}/data/genes/exons.sorted.bed.gz"
chromosome_sizes="${CONSTRAINT_TOOLS}/data/chromosome-sizes/hg19.chrom.sizes.sorted"

gaps="${CONSTRAINT_TOOLS}/data/gaps/gaps.sorted.bed.gz"
encode_exclude_regions="${CONSTRAINT_TOOLS}/data/encode-exclude-regions/encode-exclude-regions.sorted.bed.gz"

covered_sites="${CONSTRAINT_TOOLS}/data/gnomad/v3/coverage/gnomad_v3_coverage.filtered.sorted.bed.gz"

info "Compute regions putatively under neutral selection..."
info "This is the strategy Kelley Harris’s lab took" 
set -o xtrace
bedtools complement -i ${exons} -g ${chromosome_sizes} \
  | bedtools intersect -v -a - -b ${gaps} ${encode_exclude_regions} \
  | bedtools intersect -a - -b ${covered_sites} \
  | sort-compress-index-bed --name ${CONSTRAINT_TOOLS}/dist/neutral-regions
