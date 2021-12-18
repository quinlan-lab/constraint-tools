set -o errexit
set -o pipefail
set -o nounset
# set -o xtrace

source set-environment-variables.sh 

# no need to export PATH since it is already in the environment: 
# `printenv | grep -w PATH` returns non-zero output 
PATH="${CONSTRAINT_TOOLS}/install:$PATH" 

########################## 

kernel=$(uname --kernel-name)
machine=$(uname --machine)

if [[ ${machine} != 'x86_64' || ${kernel} != 'Linux'* ]]; then
  error "not Linux x86_64"
  exit 1
fi 

########################## 

if ! which conda; then 
  error "please install conda"
  exit 1
fi 

conda_environment_configuration="conda-environment.yml"

# assume that the conda environment is the first line of the conda environment configuration file
first_line=$(cat ${conda_environment_configuration} | head -1) 
IFS=: read key value <<< "${first_line}"
conda_environment_name=$(echo "${value}" | tr -d '[:space:]')

info "testing to see if a conda environment with name '${conda_environment_name}' exists..."
if conda info --envs | grep "${conda_environment_name}"; then 
  error "conda environment with name '${conda_environment_name}' already exists!"
  error "exiting..." 
  exit 1
fi 

info "creating a conda environment called ${conda_environment_name} using ${conda_environment_configuration}..."
conda env create -f ${conda_environment_configuration}

########################## 

mkdir --parents ${CONSTRAINT_TOOLS}/bin

install-jq
install-node
install-samtools
install-htslib
install-bedtools
install-bcftools

