set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset

source /scratch/ucgd/lustre-labs/quinlan/u6018199/constraint-tools/set-environment-variables.sh 

PATH="${CONSTRAINT_TOOLS}/experiments/germline-model/chen-et-al-2022:$PATH"

# https://stackoverflow.com/a/43476575/6674256
# need to export PYTHONPATH since it is not already in the environment: 
# `printenv | grep -w PYTHONPATH` returns zero output
export PYTHONPATH="${CONSTRAINT_TOOLS}/utilities"

EXTRACT_DIRECTORY="${CONSTRAINT_TOOLS_DATA}/background-selection"

download_and_preprocess () {
  rm -rf "${CONSTRAINT_TOOLS_DATA}/background-selection/CADD-B-map"

  url="https://github.com/sellalab/HumanLinkedSelectionMaps/raw/master/Bmaps/CADD_bestfit.tar.gz"
  save_path="${CONSTRAINT_TOOLS_DATA}/background-selection/CADD_bestfit.tar.gz"

  info "Downloading CADD B-map from ${url} to ${save_path}"
  wget -O "${save_path}" "${url}"

  tar -xzf "${save_path}" -C "${EXTRACT_DIRECTORY}"
  info "Extracted CADD B-map to:" "${EXTRACT_DIRECTORY}"
  rm "${save_path}"
  mv "${EXTRACT_DIRECTORY}/CADD_bestfit" "${EXTRACT_DIRECTORY}/CADD-B-map"

  # The B-maps are formatted [B, length], 
  # after McVicker et al., 2009, where B is the 
  # B value (scaled from 0-1000) covering a given 
  # chromosome segment and length is the length 
  # of the segment in base pairs. Segments sum to
  # the length of each chromosome in the hg19 
  # build.
  for chrom_id in {1..22}; do
    convert-bmap-to-bed \
      < <(\
        cat "${EXTRACT_DIRECTORY}/CADD-B-map/chr${chrom_id}.bmap.txt" \
        | awk -v id="${chrom_id}" '{print "chr"id, $1, $2}') \
      >> "${EXTRACT_DIRECTORY}/CADD-B-map/bmap.hg19.bed"
    rm "${EXTRACT_DIRECTORY}/CADD-B-map/chr${chrom_id}.bmap.txt"
    info "Converted CADD B-map for chr${chrom_id} to bed format"
  done

  sort -k1,1V -k2,2n "${EXTRACT_DIRECTORY}/CADD-B-map/bmap.hg19.bed" > "${EXTRACT_DIRECTORY}/CADD-B-map/bmap.hg19.sorted.bed"
  rm "${EXTRACT_DIRECTORY}/CADD-B-map/bmap.hg19.bed"
}

lift () {
  local filename_root=$1
  info "Lifting ${filename_root}.hg19.sorted.bed to hg38"
  bash \
    ${CONSTRAINT_TOOLS}/download-process-data/lift.sh \
    "${filename_root}.hg19.sorted.bed" \
    hg19 \
    hg38
  mv \
    ${filename_root}.hg19.sorted.bed.hg38 \
    ${filename_root}.hg38.bed
  mv \
    ${filename_root}.hg19.sorted.bed.hg38.unmapped \
    ${filename_root}.hg19.unmapped.bed
}

add_header_line () {
  echo -e "chromosome\tstart\tend\tB" \
    | cat - "${EXTRACT_DIRECTORY}/CADD-B-map/bmap.hg38.bed" \
    > "${EXTRACT_DIRECTORY}/CADD-B-map/bmap.hg38.header.bed"
}

# download_and_preprocess 
# lift "${EXTRACT_DIRECTORY}/CADD-B-map/bmap"
add_header_line





