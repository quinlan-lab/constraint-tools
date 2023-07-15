#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH --output=/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/benchmark-genome-wide-predictions/chen-et-al-2022/enhancer-characteristics-enrichment-subset-inference.log

set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset

source set-environment-variables.sh 

filename_prefix="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/enhancer-characteristics-enrichment-subset" 
windows_filename="${filename_prefix}.bed"

get-header () {
  head -1 ${windows_filename}
}

get-body () { 
  cat ${windows_filename} | tail -n +2
}

window_index=0
while read -r window; do
  echo "##############################################" 
  info "window index:" ${window_index}

  directory="${filename_prefix}-inference" 
  window_filename="${directory}/window-${window_index}.bed"

  (
    get-header
    echo "${window}"
  ) > ${window_filename}  
  
  info "Wrote:" ${window_filename}  

  info "Submitting job to infer using nucleotide transformer..."

  # %j indicates job number (JOBID) and %N indicates first node
  # use "sacct -o reqmem,maxrss,averss,elapsed -j JOBID" to determine memory usage over lifetime of job
  log_file="${directory}/window-${window_index}-%j-%N.log" 

  job_name="nucleotide-transformer.window-${window_index}"
  sbatch \
    --output=${log_file} \
    --job-name=${job_name} \
    experiments/germline-model/chen-et-al-2022/nucleotide-transformer.2.inference.window.sh \
      --directory ${directory} \
      --window-index ${window_index} 

  ((++window_index)) # increment before evaluating so that expression never evaluates to zero 
done < <(get-body)

