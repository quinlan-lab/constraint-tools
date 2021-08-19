set -o errexit
set -o pipefail
set -o nounset
set -o xtrace

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

conda_environment="constraint-tools"

instruct_to_create_conda_environment_and_activate () {
  error "please issue the following commands:" 
  info "conda create --name ${conda_environment} python=3.9" 
  info "conda activate ${conda_environment}" 
}

# https://stackoverflow.com/a/13864829/6674256
if [[ -z ${CONDA_DEFAULT_ENV+x} ]]; then 
  error "CONDA_DEFAULT_ENV is unset"
  instruct_to_create_conda_environment_and_activate
  exit 1 
else 
  info "CONDA_DEFAULT_ENV is set to '$CONDA_DEFAULT_ENV'"
fi

if [[ $CONDA_DEFAULT_ENV != ${conda_environment} ]]; then 
  error "conda environment ${conda_environment} does not exist or is not activated"
  instruct_to_create_conda_environment_and_activate
  exit 1 
fi 

pip install --requirement ${CONSTRAINT_TOOLS}/install/requirements.txt 

########################## 

mkdir --parents ${CONSTRAINT_TOOLS}/bin

install-jq
install-node
install-samtools
install-htslib
install-bedtools

