set -o errexit
set -o pipefail
# set -o noclobber
set -o xtrace
set -o nounset 

source download-data/set-environment-variables.sh 

# https://www.nature.com/articles/s41598-019-45839-z
# https://www.encodeproject.org/files/ENCFF001TDO/
encode_exclude_regions_url="https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz" 
encode_exclude_regions="encode-exclude-regions.sorted"

# data not stored at /scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/encode-exclude-regions/grch37
encode_exclude_regions_path="${CONSTRAINT_TOOLS}/data/encode-exclude-regions"

mkdir --parents ${encode_exclude_regions_path}

info "downloading and processing regions ..."
curl --location ${encode_exclude_regions_url} \
  | gunzip --stdout \
  | cut -f 1-3 \
  | get-regular-chromosomes \
  | sort-compress-index-bed --name ${encode_exclude_regions_path}/${encode_exclude_regions}
