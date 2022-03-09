set -o errexit
set -o pipefail
set -o noclobber
set -o nounset

mutations="/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/gnomad/v3/variants/gnomad_v3.sorted.tsv.gz"
genome="/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/reference/grch38/hg38.analysisSet.fa.gz"
number_chromosomes_min="140000"
kmer_size="3" 
tmpdir="/scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools/tests/germline-model/tmpdir"
window_size="51" 

mkdir --parents ${tmpdir}

source set-environment-variables.sh 

PATH="${CONSTRAINT_TOOLS}/train/germline-model:$PATH" 

# https://stackoverflow.com/a/43476575/6674256
# need to export PYTHONPATH since it is not already in the environment: 
# `printenv | grep -w PYTHONPATH` returns zero output
export PYTHONPATH="${CONSTRAINT_TOOLS}/utilities"

line_number_of_first_region="1"
line_number_of_last_region="7"
neutral_regions_for_job="${tmpdir}/neutral-regions.${line_number_of_first_region}-${line_number_of_last_region}.bed"
counts_for_job="${tmpdir}/counts.${line_number_of_first_region}-${line_number_of_last_region}.json"
# log_for_job="${tmpdir}/${line_number_of_first_region}-${line_number_of_last_region}.log"
log_for_job='stdout'

echo "" 
aggregate-counts-over-regions \
  --genome ${genome} \
  --mutations ${mutations} \
  --number-chromosomes-min ${number_chromosomes_min} \
  --kmer-size ${kmer_size} \
  --neutral-regions-filename ${neutral_regions_for_job} \
  --counts-filename ${counts_for_job} \
  --window-size ${window_size} \
  --log-filename ${log_for_job} 
