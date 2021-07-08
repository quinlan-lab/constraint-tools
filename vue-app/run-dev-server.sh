set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

export RED='\033[0;31m'
export CYAN='\033[0;36m'
export NO_COLOR='\033[0m'

root=".."

# no need to export PATH since it is already in the environment: https://unix.stackexchange.com/a/26059/406037
# this sets the path for modern versions of npm AND node
PATH="${root}/node/bin:$PATH"

bash ${root}/utilities/info.sh "npm version is: $(npm --version)"
bash ${root}/utilities/info.sh "node version is: $(node --version)"

npm run serve