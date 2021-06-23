#!/usr/bin/env bash

# https://devhints.io/bash#miscellaneous
# put option-fetching before "set -o nounset" so that we can detect flags without arguments
while [[ "$1" =~ ^- ]]; do 
  case $1 in
    --mutations ) shift; [[ ! $1 =~ ^- ]] && mutations=$1;;
    --genome ) shift; [[ ! $1 =~ ^- ]] && genome=$1;;
    --region ) shift; [[ ! $1 =~ ^- ]] && region=$1;;
    --kmer-size ) shift; [[ ! $1 =~ ^- ]] && kmer_size=$1;;
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

# create and store a list of unique tumor barcodes, if such does not exist
column_heading="Tumor_Sample_Barcode"
if [[ ! -f ${mutations%.maf.gz}.${column_heading}.txt ]]; then 
  bash ${root}/utilities/info.sh "Creating a list of unique values for the maf column: ${column_heading}"
  column_heading_index=$(python ${root}/utilities/fetch_column_heading_index.py ${mutations} ${column_heading})
  set +o errexit
  zcat ${mutations} |
    tail -n +2 | # lob off column headings
    cut -f ${column_heading_index} | # pull out column of interest 
    # head -10000 | # debug 
    sort | # required for uniq to work as expected
    uniq \
    > ${mutations%.maf.gz}.${column_heading}.txt
  set -o errexit
fi 

number_tumors=$(cat ${mutations%.maf.gz}.${column_heading}.txt | wc -l)
bash ${root}/utilities/info.sh "Number tumors: ${number_tumors}\n"

python ${root}/train-model/estimate_mutation_probabilities.py \
  --kmer-size ${kmer_size} \
  --genome ${genome} \
  --region ${region} \
  --number-tumors ${number_tumors} \
  --output ${output} \
  --mutations ${mutations} 


