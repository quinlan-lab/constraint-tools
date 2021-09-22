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

## Define directory to download files into 
var_path="${CONSTRAINT_TOOLS}/data/gnomad/v3"

## Define reference genome
reference_genome="${CONSTRAINT_TOOLS}/data/reference/grch38/hg38.analysisSet.fa.gz"

## Define python script location to process gnomad files post-download
gnomad_scripts="${CONSTRAINT_TOOLS}/download-data/gnomad/v3"

## Define chr sizes input
chr_sizes="${CONSTRAINT_TOOLS}/data/chromosome-sizes/hg38.chrom.sizes.sorted"

## Define filename for variants with sufficient coverage
coverage_file="${CONSTRAINT_TOOLS}/data/gnomad/v3/coverage/gnomad_v3_coverage.filtered.sorted.bed.gz"

mkdir --parents ${var_path}
mkdir --parents ${var_path}/intermediate_files/${chromosome}
#rm data/gnomad/v3/intermediate_files/chr1/*
#rm logs/gnomad_v3_variants/*

#######################################

## Define output names 
gnomad_variant="gnomad_v3_${chromosome}.vcf.bgz"
gnomad_vep_annotations="vep_annotations.gnomad_v3.txt"
chr_intervals="${CONSTRAINT_TOOLS}/data/gnomad/v3/intervals/hg38.chrom.intervals"

#######################################

count_running_jobs () {
        uid=$1
        prior_job_num=$2

        info "Calculating number of running processing jobs for ${uid}..."

        ## Calculate running jobs
        running_jobs=($(squeue -u ${uid} | wc -l)-${prior_job_num})
        running_jobs=$(echo ${running_jobs} | bc -l)

        info "There are currently ${running_jobs} processing jobs being run..."

        ## Return number of running jobs
        echo ${running_jobs}
}

#######################################

info "Identifying vep annotations for gnomad v3..."
bcftools view -h ${var_path}/vcf/${gnomad_variant} | grep "ID=vep" | tr ": " "\n" | grep -i Allele | sed 's/">//' | tr "|" "\n" > ${var_path}/vcf/${gnomad_vep_annotations}

info "Segmenting chromosome sizes (hg38) for processing of gnomad v3 variant file..."
chr_interval=${chr_intervals}.${chromosome}.noheader
python ${gnomad_scripts}/get_gnomad_v3_intervals.py --chr-sizes-file ${chr_sizes} --bin-num 500 --output ${chr_intervals}  
cat ${chr_intervals} | tail -n +2 | grep -w ${chromosome} | bgzip > ${chr_interval}.gz

info "Generating main.sh file to run job that extracts vcf into from variants belonging to a specific interval"
## Assign the job name
output="${CONSTRAINT_TOOLS}/logs/gnomad_v3_variants/${chromosome}/process_gnomad_v3"
job_name="${chromosome}_process_variants"

string="--gnomad-variant-file ${var_path}/vcf/${gnomad_variant} --vep-annotation-file ${var_path}/vcf/${gnomad_vep_annotations} --var-path ${var_path}"
zless ${chr_interval}.gz | awk -v gnomad_scripts="${gnomad_scripts}" -v output="${output}" -v job_name="${job_name}" -v string="${string}" '{print "sbatch --output=" output "_"$1":"$2"-"$3".out --job-name=" job_name " "  gnomad_scripts "/process_variants.sh --interval "$1":"$2"-"$3" " string}' > ${gnomad_scripts}/job_scripts/main.${chromosome}.sh
head ${gnomad_scripts}/job_scripts/main.${chromosome}.sh

#######################################

## Determine how many jobs are running prior to processing variants
info "Getting number of running jobs before all interval-specific jobs are run..."
uid="u1240855"
prior_job_num=($(squeue -u ${uid} | wc -l))

info "Executing jobs (n=many) to extract vcf info fields from variants belonging to a specific interval..."
bash ${gnomad_scripts}/job_scripts/main.${chromosome}.sh

info "Determining how many processing jobs are currently running..."
running_jobs=$(count_running_jobs ${uid} ${prior_job_num})

## Do not proceed until all processing jobs are complete
until [ ${running_jobs} -eq ${prior_job_num} ]
do
        info "Processing of gnomad v3 ${chromosome} vcf file is not complete..."
        running_jobs=$(count_running_jobs ${uid} ${prior_job_num})
done

sleep 2m

#######################################

info "Determining how many jobs were submitted..."
job_num=($(cat ${gnomad_scripts}/job_scripts/main.${chromosome}.sh | wc -l))

info "You submitted ${job_num} jobs..."

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
