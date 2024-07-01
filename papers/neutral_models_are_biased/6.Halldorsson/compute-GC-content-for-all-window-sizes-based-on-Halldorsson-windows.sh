#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-rw
#SBATCH --mem=50g # sacct -o reqmem,maxrss,averss,elapsed -j JOBID
#SBATCH --output=/scratch/ucgd/lustre-labs/quinlan/data-shared/constraint-tools/depletion_rank_scores/compute-GC-content-for-all-window-sizes-based-on-Halldorsson-windows.log

set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset

source /scratch/ucgd/lustre-labs/quinlan/u6018199/constraint-tools/set-environment-variables.sh 

# https://stackoverflow.com/a/43476575/6674256
# need to export PYTHONPATH since it is not already in the environment: 
# `printenv | grep -w PYTHONPATH` returns zero output
export PYTHONPATH="${CONSTRAINT_TOOLS}/utilities"

for window_size in "1000" "10000" "100000" "1000000"; do 
  info "Computing GC content for window size:" ${window_size}
  bash ${CONSTRAINT_TOOLS}/papers/neutral_models_are_biased/6.Halldorsson/compute-GC-content-given-window-size-based-on-Halldorsson-windows.sh ${window_size}
  info "Done computing GC content for window size:" ${window_size}
done

