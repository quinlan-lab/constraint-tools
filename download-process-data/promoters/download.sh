set -o errexit
set -o pipefail
# set -o noclobber
set -o xtrace
set -o nounset 

source download-process-data/set-environment-variables.sh 

# no need to export PATH since it is already in the environment: 
# `printenv | grep -w PATH` returns non-zero output 
PATH="${CONSTRAINT_TOOLS}/bin:$PATH" 
PATH="${CONSTRAINT_TOOLS}/utilities:$PATH" 

# https://www.sciencedirect.com/science/article/pii/S0002929720302445?via%3Dihub#app3
# Table S2
promoters_url="https://ars.els-cdn.com/content/image/1-s2.0-S0002929720302445-mmc3.csv"

promoters_path="${CONSTRAINT_TOOLS}/download-process-data/promoters"

set +o errexit
curl ${promoters_url} \
  | head -1 \
  > ${promoters_path}/promoters.headers.grch37.csv
set -o errexit

curl ${promoters_url} \
	| tail -n +2 \
  | commas-to-tabs \
	| awk -v OFS='\t' '{ print $1, $2, $3 }' \
  > ${promoters_path}/promoters.coordinates.grch37.tsv

curl ${promoters_url} \
  > ${promoters_path}/promoters.grch37.csv

