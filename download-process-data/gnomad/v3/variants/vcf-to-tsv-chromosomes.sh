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
PATH="${CONSTRAINT_TOOLS}/bin:${CONSTRAINT_TOOLS}/download-process-data/gnomad/v3/variants:$PATH"

#######################################

info "split all chromosomes into intervals ..."

chromosome_sizes="${CONSTRAINT_TOOLS_DATA}/reference/grch38/chromosome-sizes/hg38.chrom.sizes.sorted"

intervals="${CONSTRAINT_TOOLS_DATA}/map-reduce-intervals/intervals.bed"

number_of_intervals_per_job="2"
number_of_jobs="500"
# if chromosomes harbor regions of high mutation density, 
# then one may increase this number to avoid memory saturation: 
number_of_intervals_per_chromosome=$((${number_of_intervals_per_job} * ${number_of_jobs}))

split-chromosomes-into-intervals \
  --chr-sizes-file ${chromosome_sizes} \
  --number-of-intervals-per-chromosome ${number_of_intervals_per_chromosome} \
  --output ${intervals}

#######################################

# ignore X and Y to avoid complications related to allele frequencies
# TODO: uncomment
# for chromosome in $(seq 1 22); do 
for chromosome in 1; do 
  chromosome="chr${chromosome}"
  vcf-to-tsv-chromosome \
    --chromosome ${chromosome} \
    --number-of-intervals-per-job ${number_of_intervals_per_job} \
    --number-of-jobs ${number_of_jobs} \
    --number-of-intervals-per-chromosome ${number_of_intervals_per_chromosome}
done





