set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

export RED='\033[0;31m'
export CYAN='\033[0;36m'
export NO_COLOR='\033[0m'

root=$1

kernel=$(uname --kernel-name)
machine=$(uname --machine)

bash ${root}/utilities/info.sh "downloading node (and npm)"

rm -rf "${root}/node"
node_version="v16.4.2"

if [[ ${machine} == 'x86_64' ]]; then
  if [[ ${kernel} == 'Darwin'* ]]; then
    curl --remote-name "https://nodejs.org/dist/${node_version}/node-${node_version}-darwin-x64.tar.gz"
    tar -xf "node-${node_version}-darwin-x64.tar.gz"
		rm "node-${node_version}-darwin-x64.tar.gz"
		mv "node-${node_version}-darwin-x64" "${root}/node"
  elif [[ ${kernel} == 'Linux'* ]]; then
    curl --remote-name "https://nodejs.org/dist/${node_version}/node-${node_version}-linux-x64.tar.xz"
    tar -xf "node-${node_version}-linux-x64.tar.xz"
    rm "node-${node_version}-linux-x64.tar.xz"
    mv "node-${node_version}-linux-x64" "${root}/node"
  else
    echo 'neither Darwin nor Linux'
    exit 1
  fi
else
  echo 'not x86_64'
  exit 1
fi

# no need to export PATH since it is already in the environment: https://unix.stackexchange.com/a/26059/406037
# this sets the path for modern versions of npm AND node
PATH="${root}/node/bin:$PATH"

bash ${root}/utilities/info.sh "npm version is: $(npm --version)"
bash ${root}/utilities/info.sh "node version is: $(node --version)"



