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

echo -e "${CYAN}Downloading and indexing ICGC mutations...${NO_COLOR}\n"

# https://dcc.icgc.org/releases/PCAWG/consensus_snv_indel/
# Supplementary Table 4 at https://www.nature.com/articles/s41586-020-1969-6#MOESM3
url="https://dcc.icgc.org/api/v1/download?fn=/PCAWG/consensus_snv_indel/final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz"
root="/scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools"
path="${root}/data/icgc/mutations"

wget ${url} --output-document=${path}.maf.gz 
gunzip ${path}.maf.gz
bgzip ${path}.maf
tabix ${path}.maf.gz

