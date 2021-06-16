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

# p34 of: 
# https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-1969-6/MediaObjects/41586_2020_1969_MOESM1_ESM.pdf 
# says “used Ensembl release 74 as a base for the annotation”
url="http://ftp.ensembl.org/pub/release-74/gtf/homo_sapiens"
genes="Homo_sapiens.GRCh37.74"

root="/scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools"
path="${root}/data/genes"

mkdir --parents ${path}

echo -e "${CYAN}Downloading GRCH37 annotations...${NO_COLOR}\n"
wget ${url}/${genes}.gtf.gz --output-document=${path}/${genes}.gtf.gz

echo -e "${CYAN}Sorting and block compressing GRCH37 annotations...${NO_COLOR}\n"
set +o errexit
(
  zgrep "^#" ${path}/${genes}.gtf.gz
  zgrep -v "^#" ${path}/${genes}.gtf.gz | sort --version-sort -k1,1 -k4,4n
) | 
  bgzip > ${path}/${genes}.sorted.gtf.gz
set -o errexit

echo -e "${CYAN}Indexing GRCH37 annotations...${NO_COLOR}\n"
tabix ${path}/${genes}.sorted.gtf.gz
