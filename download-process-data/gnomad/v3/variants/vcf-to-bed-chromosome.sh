#!/bin/sh
#!/bin/bash
#SBATCH --time=71:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw

## Define input variables
while [[ "$1" =~ ^- ]]; do
  case $1 in
    --chromosome ) shift; [[ ! $1 =~ ^- ]] && chromosome=$1;;
  esac
  shift
done

set -o errexit
set -o pipefail
# set -o noclobber
set -o xtrace
set -o nounset 

#######################################

source download-data/set-environment-variables.sh 

#######################################

module load bcftools 
module load bgzip
module load tabix

#######################################

# https://stackoverflow.com/a/43476575/6674256
# need to export PYTHONPATH since it is not already in the environment: 
# `printenv | grep -w PYTHONPATH` returns zero output
export PYTHONPATH="${CONSTRAINT_TOOLS}/utilities:${CONSTRAINT_TOOLS}/predict-constraint"

# no need to export PATH since it is already in the environment: 
# `printenv | grep -w PATH` returns non-zero output 
PATH="${CONSTRAINT_TOOLS}/bin:$PATH"       

#######################################

## Define directory of downloaded VCFs 
var_path="${CONSTRAINT_TOOLS}/data/gnomad/v3"

## Make necessary directories
mkdir --parents ${var_path}
mkdir --parents ${var_path}/intermediate_files/${chromosome}
mkdir --parents ${var_path}/preprocess_intervals

#######################################

## Define reference genome
reference_genome="${CONSTRAINT_TOOLS}/data/reference/grch38/hg38.analysisSet.fa.gz"

## Define python script location to process gnomad files post-download
gnomad_scripts="${CONSTRAINT_TOOLS}/download-data/gnomad/v3"

## Define chr sizes input
chr_sizes="${CONSTRAINT_TOOLS}/data/chromosome-sizes/hg38.chrom.sizes.sorted"

## Define filename for variants with sufficient coverage
coverage_file="${CONSTRAINT_TOOLS}/data/gnomad/v3/coverage/gnomad_v3_coverage.filtered.sorted.bed.gz"

#######################################

## Define final processed output names 
gnomad_variant="gnomad_v3_${chromosome}.vcf.bgz"

#######################################

info "Identifying vep annotations for gnomad v3..."

## TODO: Make this its own function
gnomad_vep_annotations="${var_path}/vcf/vep_annotations.gnomad_v3.txt"
bcftools view -h ${var_path}/vcf/${gnomad_variant} | grep "ID=vep" | tr ": " "\n" | grep -i Allele | sed 's/">//' | tr "|" "\n" > ${var_path}/vcf/${gnomad_vep_annotations}

#######################################

## Define parameters for job array
job_num="500" ## This is the number of chr-specific intervals a raw vcf file will be processed on
info "Submitting job array with ${job_num} total jobs..."

#######################################

info "Defining chr-specific intervals to preprocess raw variants..."

## TODO: Make this its own function
chr_interval="${var_path}/preprocess_intervals/hg38.chrom.intervals.${chromosome}.noheader"
python ${gnomad_scripts}/get_gnomad_v3_intervals.py --chr-sizes-file ${chr_sizes} --bin-num ${job_num} --output ${chr_intervals}
cat ${chr_intervals} | tail -n +2 | grep -w ${chromosome} | bgzip > ${chr_interval}.gz

#######################################

info "Determining how many jobs were successfully completed without error..."

## Combine log output files
echo "" > ${CONSTRAINT_TOOLS}/logs/process-gnomad-v3-variants.out
for file in ${CONSTRAINT_TOOLS}/logs/gnomad_v3_variants/${chromosome}/*
do 
	tail -n 1 $file >> ${CONSTRAINT_TOOLS}/logs/process-gnomad-v3-variants.out
done

## Grep for the success signature
successful_job_num=($(cat ${CONSTRAINT_TOOLS}/logs/process-gnomad-v3-variants.out | grep "Ready to merge..." | wc -l))

if [ ${job_num} -eq ${successful_job_num} ]; then
	info "All submitted jobs to process ${chromosome} gnomad v3 vcf variants were successful... Proceeding with merge..."
	echo "" > ${var_path}/intermediate_files/gnomad_v3_variants.${chromosome}.bed
	cat ${var_path}/intermediate_files/${chromosome}/* >> ${var_path}/intermediate_files/gnomad_v3_variants.${chromosome}.bed

else 
	info "Some jobs failed... Did not merge..." 
	exit 1
fi

#######################################

info "Sorting gnomad v3 ${chromosome} variants..."
cat ${var_path}/intermediate_files/gnomad_v3_variants.${chromosome}.bed | sort -k1,1 -k2,2n > ${var_path}/intermediate_files/gnomad_v3_variants.${chromosome}.sorted.bed

info "Adding header to concatenated variant file... Also compressing..."
cat ${var_path}/header_${chromosome} | head -n 1 > ${var_path}/header
cat ${var_path}/header ${var_path}/intermediate_files/gnomad_v3_variants.${chromosome}.sorted.bed | bgzip > ${var_path}/intermediate_files/gnomad_v3_variants.${chromosome}.sorted.bed.gz

rm ${var_path}/header_${chromosome}

info "Indexing gnomad v3 variants..."
tabix \
	--skip-lines 1 \
	--sequence 1 \
	--begin 2 \
	--end 3 \
	--force \
	${var_path}/intermediate_files/gnomad_v3_variants.${chromosome}sorted.bed.gz
