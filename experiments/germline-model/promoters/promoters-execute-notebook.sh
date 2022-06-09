#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw

input_notebook="${CONSTRAINT_TOOLS}/experiments/germline-model/promoters/promoters-compute-zscores.ipynb"

number_neutral_bases="$1"

output_notebook="${CONSTRAINT_TOOLS}/experiments/germline-model/promoters/"\
"promoters-compute-zscores.${number_neutral_bases}.ipynb"

IFS=, read number_neutral_bases_lower number_neutral_bases_upper <<< "${number_neutral_bases}"

info "number_neutral_bases_lower:" ${number_neutral_bases_lower}
info "number_neutral_bases_upper:" ${number_neutral_bases_upper}

# https://papermill.readthedocs.io/en/latest/usage-inspect.html#inspect-a-notebook
# papermill --help-notebook ${input_notebook}

# https://papermill.readthedocs.io/en/latest/usage-execute.html#execute-a-notebook-with-parameters
papermill ${input_notebook} ${output_notebook} \
  -p number_neutral_bases_lower ${number_neutral_bases_lower} \
  -p number_neutral_bases_upper ${number_neutral_bases_upper}



