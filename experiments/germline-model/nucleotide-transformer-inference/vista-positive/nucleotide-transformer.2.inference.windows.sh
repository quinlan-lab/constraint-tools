#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH --output=/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/vista-enhancers/vista-enhancers.positive.hg38.hg19-inference.log

set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset

source set-environment-variables.sh 

filename_prefix="${CONSTRAINT_TOOLS_DATA}/vista-enhancers/vista-enhancers.positive.hg38.hg19" 
windows_filename="${filename_prefix}.tsv"
directory="${filename_prefix}-inference" 
mkdir -p ${directory}

get-header () {
  # head -1 ${windows_filename}
  echo -e "chromosome\tstart\tend" 
}

get-body () { 
  # cat ${windows_filename} | tail -n +2
  cat ${windows_filename} \
    | sed '/^chrX/d' \
    | sed '/^chrY/d' \
    | awk '$3 - $2 < 1000' \
    | shuf -n 100 \
    | cut -f1-3
}

window_index=0
while read -r window; do
  echo "##############################################" 
  info "window index:" ${window_index}

  window_filename="${directory}/window-${window_index}.bed"

  (
    get-header
    echo "${window}"
  ) > ${window_filename}  
  
  info "Wrote:" ${window_filename}  

  info "Submitting slurm job to infer using nucleotide transformer..."

  # %j indicates job number (JOBID) and %N indicates first node
  # use "sacct -o reqmem,maxrss,averss,elapsed -j JOBID" to determine memory usage over lifetime of job
  log_file="${directory}/window-${window_index}-%j-%N.log" 

  job_name="nucleotide-transformer.vista-positive.window-${window_index}"

  sbatch \
    --output=${log_file} \
    --job-name=${job_name} \
    experiments/germline-model/nucleotide-transformer-inference/nucleotide-transformer.2.inference.window.sh \
      --directory ${directory} \
      --window-index ${window_index} 

  ((++window_index)) # increment before evaluating so that expression never evaluates to zero 
done < <(get-body)

