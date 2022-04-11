set -o errexit
set -o pipefail
# set -o noclobber
set -o xtrace
set -o nounset 

source download-process-data/set-environment-variables.sh 

# https://academic.oup.com/bioinformatics/article/30/20/2843/2422145
# https://github.com/lh3/varcmp/blob/master/scripts/README.md
li_low_complexity_regions_url="https://github.com/lh3/varcmp/blob/master/scripts/LCR-hs38.bed.gz?raw=true"
li_low_complexity_regions="li-low-complexity-regions.sorted"

li_low_complexity_regions_path="${CONSTRAINT_TOOLS_DATA}/li-low-complexity-regions/grch38"

mkdir --parents ${li_low_complexity_regions_path}

info "downloading and processing regions ..."
curl --location ${li_low_complexity_regions_url} \
  | gunzip --stdout \
  | get-regular-chromosomes \
  | sort-compress-index-bed --name ${li_low_complexity_regions_path}/${li_low_complexity_regions}
