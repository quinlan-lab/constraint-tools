set -o errexit
set -o pipefail
set -o nounset
set -o xtrace

# https://stackoverflow.com/a/246128/6674256
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# need to export CONSTRAINT_TOOLS since it is not already in the environment: 
# `printenv | grep -w CONSTRAINT_TOOLS` returns zero output
export CONSTRAINT_TOOLS=${SCRIPT_DIR%/vue-app} # assume that the project's directory is the parent of the current-script directory

# no need to export PATH since it is already in the environment: 
# `printenv | grep -w PATH` returns non-zero output 
# this sets the path for modern versions of npm AND node
PATH="${CONSTRAINT_TOOLS}/node/bin:$PATH"
PATH="${CONSTRAINT_TOOLS}/utilities:$PATH"

info "npm version is: $(npm --version)"
info "node version is: $(node --version)"

npm run serve