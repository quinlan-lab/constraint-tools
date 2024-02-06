#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-rw
#SBATCH --mem=50g # sacct -o reqmem,maxrss,averss,elapsed -j JOBID
#SBATCH --output=experiments/germline-model/chen-et-al-2022/compute-GC-content-for-all-window-sizes-based-on-Chen-windows.log

set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset

source /scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools/set-environment-variables.sh 

PATH="${CONSTRAINT_TOOLS}/experiments/germline-model/chen-et-al-2022:$PATH"

# https://stackoverflow.com/a/43476575/6674256
# need to export PYTHONPATH since it is not already in the environment: 
# `printenv | grep -w PYTHONPATH` returns zero output
export PYTHONPATH="${CONSTRAINT_TOOLS}/utilities"

for window_size in "1000" "10000" "100000" "1000000"; do 
  info "Computing GC content for window size:" ${window_size}
  bash experiments/germline-model/chen-et-al-2022/compute-GC-content-given-window-size-based-on-Chen-windows.sh ${window_size}
  info "Done computing GC content for window size:" ${window_size}
done

