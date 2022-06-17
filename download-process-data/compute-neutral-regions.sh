set -o errexit
set -o pipefail
# set -o xtrace
set -o nounset

source download-process-data/set-environment-variables.sh 

exons="${CONSTRAINT_TOOLS_DATA}/genes/grch38/exons.sorted.bed.gz"
chromosome_sizes="${CONSTRAINT_TOOLS_DATA}/reference/grch38/chromosome-sizes/hg38.chrom.sizes.sorted"

gaps="${CONSTRAINT_TOOLS_DATA}/gaps/grch38/gaps.sorted.bed.gz"
encode_exclude_regions="${CONSTRAINT_TOOLS_DATA}/encode-exclude-regions/grch38/encode-exclude-regions.sorted.bed.gz"
test_promoters="${CONSTRAINT_TOOLS}/download-process-data/promoters/promoters.grch38.test.sorted.bed.gz" 

# giab_difficult_regions="${CONSTRAINT_TOOLS_DATA}/giab-difficult-regions/grch38/giab-difficult-regions.sorted.bed.gz"
# giab_low_complexity_regions="${CONSTRAINT_TOOLS_DATA}/giab-low-complexity-regions/grch38/giab-low-complexity-regions.sorted.bed.gz"
# li_low_complexity_regions="${CONSTRAINT_TOOLS_DATA}/li-low-complexity-regions/grch38/li-low-complexity-regions.sorted.bed.gz"

covered_sites="${CONSTRAINT_TOOLS_DATA}/gnomad/v3/gnomad_v3_coverage.filtered.sorted.bed.gz"

info "Compute regions putatively under neutral selection..."
info "This is the strategy Kelley Harris’s lab took:" 
set -o xtrace
bedtools complement -i ${exons} -g ${chromosome_sizes} \
  | bedtools intersect -v -a - -b ${gaps} ${encode_exclude_regions} ${test_promoters} \
  | bedtools intersect -a - -b ${covered_sites} \
  | sort-compress-index-bed --name ${CONSTRAINT_TOOLS}/dist/neutral-regions-germline-grch38-exclude-test-promoters
