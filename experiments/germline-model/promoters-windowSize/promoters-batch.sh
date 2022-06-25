set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset

source set-environment-variables.sh 

for window_size__window_stride in "101,10" "501,25" "1001,40" "2001,80"; do
  IFS=, read window_size window_stride <<< "${window_size__window_stride}"
  log_file="${CONSTRAINT_TOOLS}/experiments/germline-model/promoters-windowSize/promoters-batch.windowSize-${window_size}.windowStride-${window_stride}.log"
  script_file="${CONSTRAINT_TOOLS}/experiments/germline-model/promoters-windowSize/promoters-execute-notebook.sh"
  sbatch \
    --output=${log_file} \
    --job-name="promoters-batch.windowSize-${window_size}.windowStride-${window_stride}" \
    ${script_file} ${window_size__window_stride}
done