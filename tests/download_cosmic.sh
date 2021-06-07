set -o errexit
set -o pipefail
# set -o noclobber

set -o xtrace
# Must use single quote to prevent variable expansion.
# For example, if double quotes were used, ${LINENO} would take on the value of the current line,
# instead of its value when PS4 is used later in the script
# https://stackoverflow.com/a/6697845/6674256
# ${FOO:+val}    val if $FOO is set
# ${FOO[0]}   element #0 of the FOO array
# https://www.gnu.org/software/bash/manual/html_node/Bash-Variables.html
PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

export RED='\033[0;31m'
export CYAN='\033[0;36m'
export NO_COLOR='\033[0m' 

read_credentials () {
  local key_=$1
  local credentials_="cosmic_credentials"
  jq --raw-output .${key_} ${credentials_}.json
}
username=$(read_credentials "username")
password=$(read_credentials "password")
authentication_string=$(echo "${username}:${password}" | base64)

cosmic="CosmicCodingMuts.normal"
# https://cancer.sanger.ac.uk/cosmic/download?genome=38
url="https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v94/VCF/${cosmic}.vcf.gz"

authenticated_url=$(curl --header "Authorization: Basic ${authentication_string}" ${url}  | jq --raw-output .url)
curl ${authenticated_url} > ${cosmic}.vcf.gz

gunzip ${cosmic}.vcf.gz
bgzip ${cosmic}.vcf
tabix ${cosmic}.vcf.gz

