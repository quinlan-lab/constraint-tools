#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --account=redwood-gpu
#SBATCH --partition=redwood-gpu
#SBATCH --mem=20g # sacct -o reqmem,maxrss,averss,elapsed -j JOBID

# previously used: 
# --account=quinlan-rw
# --partition=quinlan-shared-rw

set -o errexit
set -o pipefail
set -o nounset
# set -o noclobber
# set -o xtrace

PATH=${CONSTRAINT_TOOLS}:$PATH 

window_size="1000"
info "Model window-size:" "${window_size}bp"
model="${CONSTRAINT_TOOLS}/dist/model-germline-grch38.windowSize-${window_size}.json"

genome_wide_predictions_directory="${CONSTRAINT_TOOLS_DATA}/genome-wide-predictions"
zscores="${genome_wide_predictions_directory}/predict-germline-grch38.exon-adjacent-windows.bed.gz"

progress_bars="disk" 
# progress_bars="stdout" 

work_directory="work-predict-germline-model-production.exon-adjacent-windows"
work_directory_should_be_clean="true"
work="${genome_wide_predictions_directory}/${work_directory}" # path to directory to store intermediate work and logs
if [[ ${work_directory_should_be_clean} == "true" && -d ${work} ]]; then 
  error "the following work directory already exists:" ${work}
  exit 1 
else 
  mkdir --parents ${work}
fi 

exon_adjacent_windows="${CONSTRAINT_TOOLS_DATA}/trustworthy-exon-adjacent-windows.bed"
filtered_exon_adjacent_windows="${CONSTRAINT_TOOLS_DATA}/filtered-trustworthy-exon-adjacent-windows.bed"

info "Windows:" ${exon_adjacent_windows}

get-filtered-exon-adjacent-windows () { 
  info "\tRemoving exon-adjacent windows on chromosomes X and Y..."
  cat ${exon_adjacent_windows} \
    | get-nonXY-chromosomes \
    > ${filtered_exon_adjacent_windows}
  info "\tWrote filtered exon-adjacent windows to:" ${filtered_exon_adjacent_windows}
}

get-filtered-exon-adjacent-windows

constraint-tools predict-germline-model \
  --model ${model} \
  --zscores ${zscores} \
  --work ${work} \
  --progress-bars ${progress_bars} \
  --number-of-jobs 500 \
  --windows ${filtered_exon_adjacent_windows}


