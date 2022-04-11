set -o errexit
set -o pipefail
# set -o noclobber
set -o xtrace
set -o nounset 

source download-process-data/set-environment-variables.sh 

# https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.0/GRCh38/LowComplexity/
# https://github.com/genome-in-a-bottle/genome-stratifications/blob/master/GRCh38/LowComplexity/GRCh38-LowComplexity-README.md
giab_low_complexity_regions_url="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.0/GRCh38/LowComplexity/GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz"
giab_low_complexity_regions="giab-low-complexity-regions.sorted"

giab_low_complexity_regions_path="${CONSTRAINT_TOOLS_DATA}/giab-low-complexity-regions/grch38"

mkdir --parents ${giab_low_complexity_regions_path}

info "downloading and processing regions ..."
curl --location ${giab_low_complexity_regions_url} \
  | gunzip --stdout \
  | get-regular-chromosomes \
  | sort-compress-index-bed --name ${giab_low_complexity_regions_path}/${giab_low_complexity_regions}
