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

reference="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/reference/GRCh38_full_analysis_set_plus_decoy_hla"

echo "***************************"
echo -e "${CYAN}test: ${NO_COLOR}"
echo -e "expected result:\nXXX"
echo -e "observed result:\nYYY"
