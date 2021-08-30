#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --nodes=1
# a slurm task is a Linux process:
# https://support.ceci-hpc.be/doc/_contents/QuickStart/SubmittingJobs/SlurmTutorial.html#going-parallel
#SBATCH --ntasks=16
# slurm does not allocate resources for more than 16 CPUs per job,
# so request one CPU per Linux process:
#SBATCH --cpus-per-task=1 
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw

set -o errexit
set -o pipefail
set -o nounset
# set -o noclobber
# set -o xtrace

CONSTRAINT_TOOLS=$1 

PATH=${CONSTRAINT_TOOLS}:$PATH 
PATH=${CONSTRAINT_TOOLS}/utilities:$PATH 
PATH=${CONSTRAINT_TOOLS}/bin:$PATH 

mutations="/scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools/data/icgc/mutations.sorted.maf.gz"
genome="/scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools/data/reference/grch37/genome.fa.gz"
kmer_size="5"
model="${CONSTRAINT_TOOLS}/dist" # path to directory to store model in

fetch_subset_of_regions () { 
  less ${CONSTRAINT_TOOLS}/dist/neutral-regions.bed.gz | head -100 
}

train_on_subset_of_regions () {
  info "$(fetch_subset_of_regions | awk '{ print $0, $3-$2 }')"

  constraint-tools train \
    --genome ${genome} \
    --mutations ${mutations} \
    --kmer-size ${kmer_size} \
    --regions <(fetch_subset_of_regions | bgzip) \
    --model ${model}
}

train_on_all_regions () {
  constraint-tools train \
    --genome ${genome} \
    --mutations ${mutations} \
    --kmer-size ${kmer_size} \
    --model ${model}
}

# train_on_subset_of_regions
train_on_all_regions