set -o errexit
set -o pipefail
# set -o noclobber
set -o xtrace
set -o nounset

source download-process-data/set-environment-variables.sh

# https://stackoverflow.com/a/43476575/6674256
# need to export PYTHONPATH since it is not already in the environment:
# `printenv | grep -w PYTHONPATH` returns zero output
export PYTHONPATH="${CONSTRAINT_TOOLS}/utilities"

# no need to export PATH since it is already in the environment:
# `printenv | grep -w PATH` returns non-zero output
PATH="${CONSTRAINT_TOOLS}/bin:$PATH"

#######################################

chromosome_sizes="${CONSTRAINT_TOOLS_DATA}/reference/grch38/chromosome-sizes/hg38.chrom.sizes.sorted"

number_of_jobs="500" 

intervals="${CONSTRAINT_TOOLS_DATA}/map-reduce-intervals/intervals.bed"

split-chromosomes-into-intervals \
  --chr-sizes-file ${chromosome_sizes} \
  --target-number-intervals ${number_of_jobs} \
  --output ${intervals}

#######################################

# TODO: loop over chromosomes and call vcf-to-bed-chromosome.sh on each 





