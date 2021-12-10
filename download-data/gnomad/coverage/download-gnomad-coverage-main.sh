#!/bin/sh
#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH --output=logs/coverage/download-gnomad-coverage.out

set -o errexit
set -o pipefail
# set -o noclobber
set -o xtrace
set -o nounset

source download-data/set-environment-variables.sh

## Define directory for log files
log_dir="${CONSTRAINT_TOOLS}/logs/coverage"

## Define file to run downloading and filtering steps
script="${CONSTRAINT_TOOLS}/download-data/gnomad/coverage/download-gnomad-coverage"

## gnomad: version 3 (WGS)
info "Downloading and filtering (by coverage) gnomad v3's WGS coverage file..."
#sbatch -W --output=${log_dir}/download-gnomad-v3-wgs-coverage.out ${script}.sh --coverage-threshold 7 --fraction-threshold 0.9 --version v3 --seq wgs

# final set of intervals in which x percent of individuals have at least y-fold coverage can be found at: 
# /scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/gnomad/v3/gnomad_v3_coverage.filtered.sorted.bed.gz


