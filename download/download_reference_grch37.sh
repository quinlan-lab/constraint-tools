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

echo -e "${CYAN}Downloading GRCH37 reference plus decoy sequences...${NO_COLOR}\n"

# https://dcc.icgc.org/releases/PCAWG/reference_data/pcawg-bwa-mem
# https://www.cureffi.org/2013/02/01/the-decoy-genome/

url="https://dcc.icgc.org/api/v1/download?fn=/PCAWG/reference_data/pcawg-bwa-mem"
root="/scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools"
path="${root}/data/reference/grch37"

mkdir --parents ${path}

download_file_with_suffix () {
  local suffix_=$1
  wget ${url}/genome.${suffix_} --output-document=${path}/genome.${suffix_}  
}

# wget ${url}/README.txt --output-document=${path}/README.txt
# download_file_with_suffix "fa.gz.fai"
# download_file_with_suffix "fa.gz"

check_digest () { 
  local suffix_=$1
  local expected_digest_=$2
  observed_digest_=$(md5sum ${path}/genome.${suffix_} | awk '{ print $1 }')
  if [[ ${observed_digest_} == ${expected_digest_} ]]; then 
    echo -e "${CYAN}${path}/genome.${suffix_}: digest check passed${NO_COLOR}"
  else
    echo -e "${RED}${path}/genome.${suffix_}: digest check failed${NO_COLOR}"
  fi 
}

check_digest "fa.gz" "a07c7647c4f2e78977068e9a4a31af15"
check_digest "fa.gz.fai" "bb77e60e9a492fd0172e2b11e6c16afd"







