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

download_basicAccessAuthentication="download_basicAccessAuthentication"
cosmic="CosmicCodingMuts.normal"
credentials="cosmic_credentials"

username=$(jq .username ${credentials}.json)
password=$(jq .password ${credentials}.json)

# auxilary file to help download COSMIC data
if [[ ! -e ${download_basicAccessAuthentication}.py ]]; then
  curl https://gist.githubusercontent.com/petermchale/f382537e680f59c169cd24c3a88c344e/raw/5dd85b2223bbe8b3c42c71ff9501c2f30ff8c5ad/${download_basicAccessAuthentication}.py > ${download_basicAccessAuthentication}.py
fi

# download COSMIC data
# https://cancer.sanger.ac.uk/cosmic/download?genome=38
echo "${CYAN}downloading, zipping, and indexing the COSMIC data file ... ${NO_COLOR}"
python ${download_basicAccessAuthentication}.py \
	--url https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v94/VCF/${cosmic}.vcf.gz \
	--username ${username} \
	--password ${password}
gunzip ${cosmic}.vcf.gz
bgzip ${cosmic}.vcf
tabix ${cosmic}.vcf.gz

