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

pip install pyyaml
conda_environment_name=$(python get_conda_environment_name.py)

if conda info --envs | grep "${conda_environment_name}"; then 
  error "conda environment with name '${conda_environment_name}' already exists!"
  exit 1
fi 

info "creating a conda environment called ${conda_environment_name} using ${conda_environment_configuration}"
conda env create -f ${conda_environment_configuration}

########################## 

mkdir --parents ${CONSTRAINT_TOOLS}/bin

install-jq
install-node
install-samtools
install-htslib
install-bedtools

