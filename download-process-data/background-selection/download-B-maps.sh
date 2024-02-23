set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset

source /scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools/set-environment-variables.sh 

PATH="${CONSTRAINT_TOOLS}/experiments/germline-model/chen-et-al-2022:$PATH"

# https://stackoverflow.com/a/43476575/6674256
# need to export PYTHONPATH since it is not already in the environment: 
# `printenv | grep -w PYTHONPATH` returns zero output
export PYTHONPATH="${CONSTRAINT_TOOLS}/utilities"

rm -rf "${CONSTRAINT_TOOLS_DATA}/background-selection/CADD-B-map"

url="https://github.com/sellalab/HumanLinkedSelectionMaps/raw/master/Bmaps/CADD_bestfit.tar.gz"
save_path="${CONSTRAINT_TOOLS_DATA}/background-selection/CADD_bestfit.tar.gz"

info "Downloading CADD B-map from ${url} to ${save_path}"
wget -O "${save_path}" "${url}"

extract_directory="${CONSTRAINT_TOOLS_DATA}/background-selection"
tar -xzf "${save_path}" -C "${extract_directory}"
info "Extracted CADD B-map to:" "${extract_directory}"
rm "${save_path}"
mv "${extract_directory}/CADD_bestfit" "${extract_directory}/CADD-B-map"

# The B-maps are formatted [B, length], 
# after McVicker et al., 2009, where B is the 
# B value (scaled from 0-1000) covering a given 
# chromosome segment and length is the length 
# of the segment in base pairs. Segments sum to
# the length of each chromosome in the hg19 
# build.

# convert B-maps to bed format: 
for chrom_id in {1..22}; do
  convert-bmap-to-bed \
    < <(\
      cat "${extract_directory}/CADD-B-map/chr${chrom_id}.bmap.txt" \
      | awk -v id="${chrom_id}" '{print "chr"id, $1, $2}') \
    >> "${extract_directory}/CADD-B-map/bmap.hg19.bed"
  rm "${extract_directory}/CADD-B-map/chr${chrom_id}.bmap.txt"
  info "Converted CADD B-map for chr${chrom_id} to bed format"
done

