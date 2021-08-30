set -o errexit
set -o pipefail
# set -o noclobber

set -o xtrace

source download-data/set-environment-variables.sh 

#chromosome_sizes_url="ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips"
#chromosome_sizes="hg19.chrom.sizes"

chromosome_sizes_url="https://hgdownload-test.gi.ucsc.edu/goldenPath/hg38/bigZips"
chromosome_sizes="hg38.chrom.sizes"

chromosome_sizes_path="${CONSTRAINT_TOOLS}/data/chromosome-sizes"

mkdir --parents ${chromosome_sizes_path}

info "Downloading chromosome sizes..."
curl --location ${chromosome_sizes_url}/${chromosome_sizes} \
  | get-regular-chromosomes \
  | sort --version-sort -k1,1 -k2,2 \
  > ${chromosome_sizes_path}/${chromosome_sizes}.sorted
