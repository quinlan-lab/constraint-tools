set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset

source download-process-data/set-environment-variables.sh 

# no need to export PATH since it is already in the environment: 
# `printenv | grep -w PATH` returns non-zero output 
PATH="${CONSTRAINT_TOOLS}/bin:$PATH" 
PATH="${CONSTRAINT_TOOLS}/utilities:$PATH" 

# exporting these variables make them visible in python scripts below 
export CDTS_PATH="${CONSTRAINT_TOOLS_DATA}/CDTS"

info "CDTS directory:" ${CDTS_PATH}

uncompress_and_reformat () { 
  bgzip -d -c ${CDTS_PATH}/CDTS_diff_perc_coordsorted_gnomAD_N15496_hg19.bed.gz \
    | process_CDTS_data \
    > ${CDTS_PATH}/CDTS.gnomAD.hg19.bed
}

lift () {
  info "Lifting over..."
  bash \
    ${CONSTRAINT_TOOLS}/download-process-data/lift.sh \
    ${CDTS_PATH}/CDTS.gnomAD.hg19.bed \
    hg19 \
    hg38 \
    "-bedPlus=3 -tab"
  mv \
    "${CDTS_PATH}/CDTS.gnomAD.hg19.bed.hg38" \
    "${CDTS_PATH}/CDTS.gnomAD.hg38.bed"
  mv \
    "${CDTS_PATH}/CDTS.gnomAD.hg19.bed.hg38.unmapped" \
    "${CDTS_PATH}/CDTS.gnomAD.hg19.umapped.bed"
}

print_header_line_for_CDTS () {
  # /scratch/ucgd/lustre-labs/quinlan/data-shared/constraint-tools/CDTS/README
  echo -e "chrom\tstart\tend\tobserved_counts\texpected_counts\tobserved_minus_expected\tpercentile_rank_of_observed_minus_expected"
}

add_header_line () {
  print_header_line_for_CDTS \
    | cat - "${CDTS_PATH}/CDTS.gnomAD.hg38.bed" \
    > "${CDTS_PATH}/CDTS.gnomAD.hg38.header.bed"
}

# uncompress_and_reformat
# lift 
add_header_line
