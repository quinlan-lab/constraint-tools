#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH --output=download-process-data/download-constitutive-exons-grch38.log

set -o errexit
set -o pipefail
# set -o noclobber

set -o xtrace
set -o nounset 

source download-process-data/set-environment-variables.sh 

# https://stackoverflow.com/a/43476575/6674256
# need to export PYTHONPATH since it is not already in the environment: 
# `printenv | grep -w PYTHONPATH` returns zero output
export PYTHONPATH="${CONSTRAINT_TOOLS}/utilities"

script_path=${CONSTRAINT_TOOLS}/download-process-data

genes_path="${CONSTRAINT_TOOLS_DATA}/genes/grch38"
mkdir --parents ${genes_path}

python ${script_path}/download-biomart.py ${script_path}/constitutive-exons.xml \
  | get-regular-chromosomes \
  | awk --assign OFS='\t' '$5 == "1" { print "chr"$1, $2, $3, $4 }' \
  | sort-compress-index-bed --name ${genes_path}/constitutive-exons.sorted
