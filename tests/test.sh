set -o errexit
set -o pipefail
# set -o noclobber

# set -o xtrace
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

data_root="/scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools"
mutations="${data_root}/data/icgc/mutations.sorted.maf.gz"
genome="${data_root}/data/reference/grch37/genome.fa.gz"

neutral_region="22:30,000,000-31,000,000"
kmer_size="5"

output="$PWD/tests" # assumes that this script is run from user's constraint-tools directory

PATH="$PWD:$PATH" # assumes that this script is run from user's constraint-tools directory

constraint-tools train \
  --genome ${genome} \
  --region ${neutral_region} \
  --mutations ${mutations} \
  --kmer-size ${kmer_size} \
  --output ${output}

model="${output}/model.json"
port="5000"

constraint-tools visualize \
  --model ${model} \
  --port ${port}
