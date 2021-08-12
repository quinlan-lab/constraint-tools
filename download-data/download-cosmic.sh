set -o errexit
set -o pipefail
# set -o noclobber
set -o xtrace 

source set-environment-variables.sh 

read_credentials () {
  local key_=$1
  local credentials_="cosmic-credentials"
  jq --raw-output .${key_} ${credentials_}.json
}
username=$(read_credentials "username")
password=$(read_credentials "password")
authentication_string=$(echo "${username}:${password}" | base64)

cosmic="CosmicCodingMuts.normal"
# https://cancer.sanger.ac.uk/cosmic/download?genome=38
cosmic_url="https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v94/VCF/${cosmic}.vcf.gz"

authenticated_url=$(curl --header "Authorization: Basic ${authentication_string}" ${cosmic_url}  | jq --raw-output .url)
info "Downloading COSMIC...\n"
cosmic_path="${CONSTRAINT_TOOLS}/data/cosmic/${cosmic}"
curl ${authenticated_url} > ${cosmic_path}.vcf.gz 

gunzip ${cosmic_path}.vcf.gz
bgzip ${cosmic_path}.vcf
tabix ${cosmic_path}.vcf.gz

