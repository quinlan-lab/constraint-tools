#!/usr/bin/env bash

# https://devhints.io/bash#miscellaneous
# put option-fetching before "set -o nounset" so that we can detect flags without arguments
while [[ "$1" =~ ^- ]]; do 
  case $1 in
    --kmer-size ) shift; [[ ! $1 =~ ^- ]] && kmer_size=$1;;
    --chromosome ) shift; [[ ! $1 =~ ^- ]] && chromosome=$1;;
    --start ) shift; [[ ! $1 =~ ^- ]] && start=$1;;
    --end ) shift; [[ ! $1 =~ ^- ]] && end=$1;;
    --genome ) shift; [[ ! $1 =~ ^- ]] && genome=$1;;
    --model ) shift; [[ ! $1 =~ ^- ]] && model=$1;;
    --number-tumors ) shift; [[ ! $1 =~ ^- ]] && number_tumors=$1;;
    --mutations ) shift; [[ ! $1 =~ ^- ]] && mutations=$1;;
    --output ) shift; [[ ! $1 =~ ^- ]] && output=$1;;
    --root ) shift; [[ ! $1 =~ ^- ]] && root=$1;;
    *) echo -e "${RED}$0: $1 is an invalid flag${NO_COLOR}" >&2; exit 1;;
  esac 
  shift
done

set -o errexit
set -o pipefail
set -o noclobber
set -o nounset

# set -o xtrace
# Must use single quote to prevent variable expansion.
# For example, if double quotes were used, ${LINENO} would take on the value of the current line,
# instead of its value when PS4 is used later in the script
# https://stackoverflow.com/a/6697845/6674256
# ${FOO:+val}    val if $FOO is set
# ${FOO[0]}   element #0 of the FOO array
# https://www.gnu.org/software/bash/manual/html_node/Bash-Variables.html
PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

## Count the number of tumors available as done in train-model.sh
number_tumors=$(cat ${mutations%.maf.gz}.${column_heading}.txt | wc -l)
bash ${root}/utilities/info.sh "Number tumors: ${number_tumors}\n"
  
python ${root}/v2_scripts/calculate_purifying_probabilities.py \
  --kmer-size ${kmer_size} \
  --chromosome ${chromosome} \
  --start ${start} \
  --end ${end} \
  --genome ${genome} \
  --model ${model} \
  --number-tumors ${number-tumors} \
  --mutations ${mutations} \
  --output ${output}


