set -o errexit
set -o pipefail
# set -o noclobber
set -o xtrace
set -o nounset 

source set-environment-variables.sh 

# no need to export PATH since it is already in the environment: 
# `printenv | grep -w PATH` returns non-zero output 
PATH="${CONSTRAINT_TOOLS}/experiments/germline-model/promoters-windowSize:$PATH" 

train-test-split-promoters-core

stream-test-promoter-coordinates \
  | sort-compress-index-bed --name ${CONSTRAINT_TOOLS}/download-process-data/promoters/promoters.grch38.test.sorted
