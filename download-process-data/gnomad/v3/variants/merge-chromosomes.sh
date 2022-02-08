#!/usr/bin/env bash
#SBATCH --time=4:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH --output=/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/gnomad/v3/variants/merge-chromosomes.log

set -o errexit
set -o pipefail
# set -o noclobber
set -o xtrace
set -o nounset

source download-process-data/set-environment-variables.sh

# no need to export PATH since it is already in the environment:
# `printenv | grep -w PATH` returns non-zero output
PATH="${CONSTRAINT_TOOLS}/bin:${CONSTRAINT_TOOLS}/utilities:$PATH"

variants_directory="${CONSTRAINT_TOOLS_DATA}/gnomad/v3/variants" 

random_chromosome_file=$(/usr/bin/ls ${variants_directory}/*.tsv.gz | head -1) || true

merge-chromosomes () {
  for i in $(seq 1 22); do
    set +o errexit 
    zcat ${variants_directory}/gnomad_v3_chr${i}.sorted.tsv.gz | tail -n +2 
    set -o errexit 
  done 
} 

info "Merging tsv files for individual chromosomes, and block compressing, ..."
(
  set +o errexit
  zcat ${random_chromosome_file} | head -1
  set -o errexit 
  merge-chromosomes
) | bgzip > ${variants_directory}/gnomad_v3.sorted.tsv.gz

info "Indexing tsv file for all chromosomes..."
# http://www.htslib.org/doc/tabix.html
tabix \
    --skip-lines 1 \
    --sequence 1 \
    --begin 2 \
    --end 3 \
    --force \
  ${variants_directory}/gnomad_v3.sorted.tsv.gz

