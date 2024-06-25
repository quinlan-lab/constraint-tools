set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset

source download-process-data/set-environment-variables.sh 

# no need to export PATH since it is already in the environment: 
# `printenv | grep -w PATH` returns non-zero output 
PATH="${CONSTRAINT_TOOLS}/bin:$PATH" 
PATH="${CONSTRAINT_TOOLS}/utilities:$PATH" 

# https://www.nature.com/articles/s41586-022-04965-x
# GRCH38

# info "Downloading..." 
# curl https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-04965-x/MediaObjects/41586_2022_4965_MOESM3_ESM.gz \
#   > ${CONSTRAINT_TOOLS_DATA}/depletion_rank_scores/41586_2022_4965_MOESM3_ESM.gz
# info "Finished downloading." 

gunzip ${CONSTRAINT_TOOLS_DATA}/depletion_rank_scores/41586_2022_4965_MOESM3_ESM.gz