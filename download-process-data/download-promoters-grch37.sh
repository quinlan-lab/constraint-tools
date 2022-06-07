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

promoters_url="https://ars.els-cdn.com/content/image/1-s2.0-S0002929720302445-mmc3.csv"
promoters_path="${CONSTRAINT_TOOLS_DATA}/promoters/grch37"
mkdir --parents ${promoters_path}

# https://www.sciencedirect.com/science/article/pii/S0002929720302445?via%3Dihub#app3
# Table S2
curl ${promoters_url} > ${promoters_path}/promoters.csv

