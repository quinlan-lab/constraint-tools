
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

# gap file: ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz
# gap file contains telomere and centromere information


url="https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz"
root="/scratch/ucgd/lustre-work/quinlan/u1240855/constraint-tools"
path="${root}/data/encode"

mkdir --parents ${path}

download_file () {
    local url=$1
    wget ${url} --output-document=${path}/encode_blacklist.bed.gz
}

echo -e "${CYAN}Downloading UCSC gap file...${NO_COLOR}\n"
download_file ${url}

echo -e "${CYAN}Decompressing UCSC gap file...${NO_COLOR}\n"
set +o errexit
gzip --decompress ${path}/encode_blacklist.bed.gz
set -o errexit
