set -o errexit
set -o pipefail
set -o nounset
# set -o xtrace
PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

export RED='\033[0;31m'
export CYAN='\033[0;36m'
export NO_COLOR='\033[0m'

# https://stackoverflow.com/a/246128/6674256
root="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

########################## 

kernel=$(uname --kernel-name)
machine=$(uname --machine)

if [[ ${machine} != 'x86_64' || ${kernel} != 'Linux'* ]]; then
  bash ${root}/utilities/error.sh "not Linux x86_64"
  exit 1
fi 

########################## 

if ! which conda; then 
  bash ${root}/utilities/error.sh "please install conda"
  exit 1
fi 

conda_environment="constraint-tools"
if [[ $CONDA_DEFAULT_ENV != ${conda_environment} ]]; then 
  bash ${root}/utilities/error.sh "conda environment ${conda_environment} does not exist or is not activated"
  bash ${root}/utilities/error.sh "please issue the following commands:" 
  bash ${root}/utilities/info.sh "conda create --name ${conda_environment} python=3.9" 
  bash ${root}/utilities/info.sh "conda activate ${conda_environment}" 
  exit 1 
fi 

pip install --requirement ${root}/install/requirements.txt 

########################## 

mkdir --parents ${root}/bin

bash install/jq.sh ${root}

########################## 

bash install/node.sh ${root}


