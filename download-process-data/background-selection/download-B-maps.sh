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

url="https://github.com/sellalab/HumanLinkedSelectionMaps/raw/master/Bmaps/CADD_bestfit.tar.gz"
save_path="${CONSTRAINT_TOOLS_DATA}/background-selection/CADD_bestfit.tar.gz"

info "Downloading CADD bestfit data from ${url} to ${save_path}"
wget -O "${save_path}" "${url}"

# Extract and uncompress the .tar.gz file
extract_directory="${CONSTRAINT_TOOLS_DATA}/background-selection"
tar -xzf "${save_path}" -C "${extract_directory}"
info "Extracted CADD bestfit data to:" "${extract_directory}"
rm "${save_path}"

