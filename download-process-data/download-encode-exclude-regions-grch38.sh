set -o errexit
set -o pipefail
# set -o noclobber
set -o xtrace
set -o nounset 

source download-data/set-environment-variables.sh 

# https://www.nature.com/articles/s41598-019-45839-z
# https://www.encodeproject.org/files/ENCFF356LFX/
encode_exclude_regions_url="https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz" 
encode_exclude_regions="encode-exclude-regions.sorted"

encode_exclude_regions_path="${CONSTRAINT_TOOLS_DATA}/encode-exclude-regions/grch38"

mkdir --parents ${encode_exclude_regions_path}

info "downloading and processing regions ..."
curl --location ${encode_exclude_regions_url} \
  | gunzip --stdout \
  | get-regular-chromosomes \
  | sort-compress-index-bed --name ${encode_exclude_regions_path}/${encode_exclude_regions}
