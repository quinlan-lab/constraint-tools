set -o errexit
set -o pipefail
# set -o noclobber
set -o xtrace
set -o nounset 

source download-process-data/set-environment-variables.sh 

# https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.0/GRCh38/union/GRCh38-union-README.md
giab_difficult_regions_url="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.0/GRCh38/union/GRCh38_alldifficultregions.bed.gz"
giab_difficult_regions="giab-difficult-regions.sorted"

giab_difficult_regions_path="${CONSTRAINT_TOOLS_DATA}/giab-difficult-regions/grch38"

mkdir --parents ${giab_difficult_regions_path}

info "downloading and processing regions ..."
curl --location ${giab_difficult_regions_url} \
  | gunzip --stdout \
  | get-regular-chromosomes \
  | sort-compress-index-bed --name ${giab_difficult_regions_path}/${giab_difficult_regions}
