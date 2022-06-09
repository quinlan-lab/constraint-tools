set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset

source set-environment-variables.sh 

for number_neutral_bases in "150,200" "200,250" "250,300" "300,400" "400,500" "500,600" "600,750" "750,1000"; do
  log_file="${CONSTRAINT_TOOLS}/experiments/germline-model/promoters/promoters-batch.${number_neutral_bases}.log"
  script_file="${CONSTRAINT_TOOLS}/experiments/germline-model/promoters/promoters-execute-notebook.sh"
  sbatch \
    --output=${log_file} \
    --job-name="promoters-batch.${number_neutral_bases}" \
    ${script_file} ${number_neutral_bases}
done