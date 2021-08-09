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


root="/scratch/ucgd/lustre-work/quinlan/u1240855/constraint-tools"
genome_file="${root}/data/reference/grch37/genome.fa.gz"
mutation_file="${root}/data/icgc/mutations.sorted.maf.gz"
neutral_regions="${root}/data/neutral_region/putative_neutral_regions_noprefix_numchr.bed"
model="${root}/v2_scripts/output/model.json"
kmer_size="3"
output="/scratch/ucgd/lustre-work/quinlan/u1240855/constraint-tools/v2_scripts/output"
chromosome="14"
start="102084858"
end="102084905"

echo ${root}

#bash ${root}/v2_scripts/train-model.sh --mutations ${mutation_file} --genome ${genome_file} --neutral_regions ${neutral_regions} --kmer-size ${kmer_size} --output ${output} --root ${root}

bash ${root}/v2_scripts/predict-model.sh --kmer-size ${kmer_size} --chromosome ${chromosome} --start ${start} --end ${end} \
--genome ${genome_file} --model ${model} --mutations ${mutation_file} --output ${output} --root ${root}




