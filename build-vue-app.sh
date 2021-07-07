#!/usr/bin/env bash

set -o errexit
set -o pipefail
set -o noclobber
set -o nounset

set -o xtrace
# Must use single quote to prevent variable expansion.
# For example, if double quotes were used, ${LINENO} would take on the value of the current line,
# instead of its value when PS4 is used later in the script
# https://stackoverflow.com/a/6697845/6674256
# ${FOO:+val}    val if $FOO is set
# ${FOO[0]}   element #0 of the FOO array
# https://www.gnu.org/software/bash/manual/html_node/Bash-Variables.html
PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

root=$PWD
cd ${root}/vue-app
npm install
npm run build 
cd ${root}
rm --recursive ${root}/flask-app/static
mkdir --parents ${root}/flask-app/static
mv ${root}/vue-app/dist/* ${root}/flask-app/static