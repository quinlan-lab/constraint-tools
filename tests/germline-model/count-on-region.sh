set -o errexit
set -o pipefail
set -o noclobber
set -o nounset

mutations="/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/gnomad/v3/variants/gnomad_v3.sorted.tsv.gz"
genome="/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/reference/grch38/hg38.analysisSet.fa.gz"
neutral_region="chr1:15363-15768"
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

count-on-region \
  --genome ${genome} \
  --mutations ${mutations} \
  --number-chromosomes-min ${number_chromosomes_min} \
  --kmer-size ${kmer_size} \
  --tmpdir ${tmpdir} \
  --neutral-region ${neutral_region} \
  --window-size ${window_size}

