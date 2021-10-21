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

## Run sbatch scripts to download and filter gnomad coverage files

## gnomad: version 3 (WGS)
info "Downloading and filtering gnomad v3's WGS coverage file..."
#sbatch -W --output=${log_dir}/download-gnomad-v3-wgs-coverage.out ${script}.sh --coverage-threshold 7 --fraction-threshold 0.9 --version v3 --seq wgs

## gnomad: version 2 WGS
info "Downloading and fitering gnomad v2's WGS coverage file..."
#sbatch -W --output=${log_dir}/download-gnomad-v2-wgs-coverage.out ${script}.sh --coverage-threshold 7 --fraction-threshold 0.9 --version v2 --seq wgs

## gnomad: version 2 WES
info "DOwnloading and filtering gnomad v2's WES coverage file..."
#sbatch -W --output=${log_dir}/download-gnomad-v2-wes-coverage.out ${script}.sh --coverage-threshold 7 --fraction-threshold 0.9 --version v2 --seq wes

#######################################

info "Combining gnomad v2 and v3 coverage files..."
## Define concatenated coverage file
combined_wgs_cov_file=${CONSTRAINT_TOOLS}/dist/covered_sites_wgs.bed
combined_wes_cov_file=${CONSTRAINT_TOOLS}/dist/covered_sites_wes.bed


## Define coverage files to combine
gnomad_v3="${CONSTRAINT_TOOLS}/data/gnomad/v3/coverage/gnomad_v3_coverage.filtered.sorted.bed.gz"
gnomad_v2_wgs="${CONSTRAINT_TOOLS}/data/gnomad/v2/wgs/coverage/gnomad_v2_wgs_coverage.filtered.sorted.bed.gz"
gnomad_v2_wes="${CONSTRAINT_TOOLS}/data/gnomad/v2/wes/coverage/gnomad_v2_wes_coverage.filtered.sorted.bed.gz"

## Run bedtools intersect on the above files to identify regions covered in all datasets
info "Generating combined coverage file in exonic space..."
bedtools intersect -sorted -a ${gnomad_v2_wes} -b ${gnomad_v3} ${gnomad_v2_wgs} | sort -k1,1 -k2,2n > ${combined_wes_cov_file}

info "Generating combined coverage file across the entire genome (including intergenic space)..."
bedtools intersect -sorted -a ${gnomad_v3} -b ${gnomad_v2_wgs} | sort -k1,1 -k2,2n > ${combined_wgs_cov_file}



