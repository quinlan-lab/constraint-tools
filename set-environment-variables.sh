set -o errexit 

# https://stackoverflow.com/a/246128/6674256
# need to export CONSTRAINT_TOOLS since it is not already in the environment: 
# `printenv | grep -w CONSTRAINT_TOOLS` returns zero output
export CONSTRAINT_TOOLS="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# no need to export PATH since it is already in the environment: 
# `printenv | grep -w PATH` returns non-zero output 
PATH="${CONSTRAINT_TOOLS}/utilities:$PATH" 

# info "${CONSTRAINT_TOOLS}/set-environment-variables.sh: $(printenv | grep CONSTRAINT_TOOLS)"
# info "${CONSTRAINT_TOOLS}/set-environment-variables.sh: $PATH"
