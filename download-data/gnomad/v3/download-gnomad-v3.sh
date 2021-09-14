#!/bin/sh
#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH -o logs/download-gnomad-v3.out

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

## Define chromosome to download
#chromosome="chr1"

#######################################

## Define files to download
gnomad_url="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.1/vcf/genomes/gnomad.genomes.v3.1.1.sites.${chromosome}.vcf.bgz"
gnomad_tbi_url="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.1/vcf/genomes/gnomad.genomes.v3.1.1.sites.${chromosome}.vcf.bgz.tbi"

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
#rm download-data/gnomad/v3/job_scripts/main.chr1.sh
#rm data/gnomad/v3/intervals/hg38.chrom.intervals.chr1.noheader.gz

#######################################

## Define output names 
gnomad_variant="gnomad_v3_${chromosome}.vcf.bgz"
gnomad_vep_annotations="vep_annotations.gnomad_v3.txt"
chr_intervals="${CONSTRAINT_TOOLS}/data/gnomad/v3/intervals/hg38.chrom.intervals"

#######################################

info "Downloading gnomad's variant VCF file..."
wget ${gnomad_url} --output-document=${var_path}/vcf/${gnomad_variant}

info "Downloading gnomad's tbi file..."
wget ${gnomad_tbi_url} --output-document=${var_path}/vcf/${gnomad_variant}.tbi

#######################################

info "Identifying vep annotations for gnomad v3..."
#bcftools view -h ${var_path}/vcf/${gnomad_variant} | grep "ID=vep" | tr ": " "\n" | grep -i Allele | sed 's/">//' | tr "|" "\n" > ${var_path}/vcf/${gnomad_vep_annotations}

info "Segmenting chromosome sizes (hg38) for processing of gnomad v3 variant file..."
chr_interval=${chr_intervals}.${chromosome}.noheader
python ${gnomad_scripts}/get_gnomad_v3_intervals.py --chr-sizes-file ${chr_sizes} --bin-num 500 --output ${chr_intervals}  
cat ${chr_intervals} | tail -n +2 | grep -w ${chromosome} | bgzip > ${chr_interval}.gz

info "Generating main.sh file to run job that extracts vcf into from variants belonging to a specific interval"
string="--chromosome ${chromosome} --gnomad-variant-file ${var_path}/vcf/${gnomad_variant} --vep-annotation-file ${var_path}/vcf/${gnomad_vep_annotations} --var-path ${var_path}"
zless ${chr_interval}.gz | awk -v gnomad_scripts="${gnomad_scripts}" -v string="${string}" '{print "sbatch " gnomad_scripts "/process_variants.sh --interval "$1":"$2"-"$3" " string}' > ${gnomad_scripts}/job_scripts/main.${chromosome}.sh

info "Executing jobs (n=many) to extract vcf info fields from variants belonging to a specific interval..."
bash ${gnomad_scripts}/job_scripts/main.${chromosome}.sh

info "Waiting for all jobs to complete before proceeding..."
wait

info "Determining how many jobs were submitted..."
job_num=($(cat ${gnomad_scripts}/job_scripts/main.${chromosome}.sh | wc -l))

info "You submitted ${job_num} jobs..."

info "Determining how many jobs were successfully completed without error..."
## Combine log output files
echo "" > ${CONSTRAINT_TOOLS}/logs/process-gnomad-v3-variants.out
for file in ${CONSTRAINT_TOOLS}/logs/gnomad_v3_variants/*
do 
	tail -n 1 $file >> ${CONSTRAINT_TOOLS}/logs/process-gnomad-v3-variants.out
done

## Grep for the success signature
successful_job_num=($(cat ${CONSTRAINT_TOOLS}/logs/process-gnomad-v3-variants.out | grep "Ready to merge..." | wc -l))

if [ ${job_num} -eq ${successful_job_num} ]; then
	info "All submitted jobs to process ${chromosome} gnomad v3 vcf variants were successful... Proceeding with merge..."
	echo "" > ${var_path}/intermediate_files/gnomad_v3_variants.${chromosome}.bed
	cat ${var_path}/intermediate_files/${chromosome}/* >> cat ${var_path}/intermediate_files/gnomad_v3_variants.${chromosome}.bed

else 
	info "Some jobs failed... Did not merge..." 
	exit 1
fi

info "Adding header to concatenated variant file..."
cat ${var_path}/header_${chromosome} | head -n 1 > ${var_path}/header_${chromosome}
cat ${var_path}/header_${chromosome} ${var_path}/intermediate_files/gnomad_v3_variants.${chromosome}.bed > ${var_path}/intermediate_files/gnomad_v3_variants.${chromosome}.header.bed

info "Sorting and compressing gnomad v3 variants..."
cat ${var_path}/intermediate_files/gnomad_v3_variants.${chromosome}.header.bed | tail -n +2 | sort -k1,1 -k2,2n | bgzip > ${var_path}/intermediate_files/gnomad_v3_variants.${chromosome}.sorted.bed.gz

rm ${var_path}/intermediate_files/gnomad_v3_variants.${chromosome}.header.bed

info "Indexing gnomad v3 variants..."
tabix \
	--skip-lines 1 \
	--sequence 1 \
	--begin 2 \
	--end 3 \
	--force \
	${var_path}/intermediate_files/gnomad_v3_variants.${chromosome}sorted.bed.gz
