#!/usr/bin/env bash

set -o errexit
set -o pipefail
set -o noclobber
set -o nounset

set -o xtrace

source set-environment-variables.sh 

# This sets the path to modern versions of npm AND node.
# No need to export PATH since it is already in the environment: 
# `printenv | grep -w PATH` returns non-zero output 
PATH="${CONSTRAINT_TOOLS}/node/bin:$PATH" 

info "npm version is: $(npm --version)"
info "node version is: $(node --version)"

cd ${CONSTRAINT_TOOLS}/vue-app
npm install
npm run build 
cd ${CONSTRAINT_TOOLS}
static=${CONSTRAINT_TOOLS}/flask-app/germline-model/static
rm --force --recursive ${static} 
mkdir --parents ${static} 
mv ${CONSTRAINT_TOOLS}/vue-app/dist/* ${static} 
