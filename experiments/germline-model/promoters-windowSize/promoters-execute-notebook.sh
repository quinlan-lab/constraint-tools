#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw

input_notebook="${CONSTRAINT_TOOLS}/experiments/germline-model/promoters-windowSize/compute-zscores-on-test-promoters.ipynb"

window_size__window_stride="$1"

IFS=, read window_size window_stride <<< "${window_size__window_stride}"

info "window_size:" ${window_size}
info "window_stride:" ${window_stride}

output_notebook="${CONSTRAINT_TOOLS}/experiments/germline-model/promoters-windowSize/"\
"compute-zscores-on-test-promoters.windowSize-${window_size}.windowStride-${window_stride}.ipynb"

# https://papermill.readthedocs.io/en/latest/usage-inspect.html#inspect-a-notebook
# papermill --help-notebook ${input_notebook}

# https://papermill.readthedocs.io/en/latest/usage-execute.html#execute-a-notebook-with-parameters
papermill ${input_notebook} ${output_notebook} \
  -p window_size ${window_size} \
  -p window_stride ${window_stride}



