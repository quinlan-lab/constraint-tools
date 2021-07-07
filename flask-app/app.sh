#!/usr/bin/env bash

# https://devhints.io/bash#miscellaneous
# put option-fetching before "set -o nounset" so that we can detect flags without arguments
while [[ "$1" =~ ^- ]]; do 
  case $1 in
    --model ) shift; [[ ! $1 =~ ^- ]] && model=$1;;
    --port ) shift; [[ ! $1 =~ ^- ]] && port=$1;;
    --root ) shift; [[ ! $1 =~ ^- ]] && root=$1;;
    *) echo -e "${RED}$0: $1 is an invalid flag${NO_COLOR}" >&2; exit 1;;
  esac 
  shift
done

set -o errexit
set -o pipefail
set -o noclobber
set -o nounset

# Must use single quote to prevent variable expansion.
# For example, if double quotes were used, ${LINENO} would take on the value of the current line,
# instead of its value when PS4 is used later in the script
# https://stackoverflow.com/a/6697845/6674256
# ${FOO:+val}    val if $FOO is set
# ${FOO[0]}   element #0 of the FOO array
# https://www.gnu.org/software/bash/manual/html_node/Bash-Variables.html
PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

cd ${root}/flask-app
FLASK_ENV="development" # https://flask.palletsprojects.com/en/2.0.x/server/#in-code
python app.py \
  --model ${model} \
  --port ${port}


