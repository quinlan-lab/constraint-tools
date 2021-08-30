#!/bin/sh
#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --nodes=3
#SBATCH --ntasks=16
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH -o logs/download-gnomad-v3.out

set -o errexit
set -o pipefail
# set -o noclobber
set -o xtrace
set -o nounset 

#######################################

source download-data/set-environment-variables.sh 

#######################################

module load bcftools 

#######################################

# https://stackoverflow.com/a/43476575/6674256
# need to export PYTHONPATH since it is not already in the environment: 
# `printenv | grep -w PYTHONPATH` returns zero output
export PYTHONPATH="${CONSTRAINT_TOOLS}/utilities:${CONSTRAINT_TOOLS}/predict-constraint"

# no need to export PATH since it is already in the environment: 
# `printenv | grep -w PATH` returns non-zero output 
PATH="${CONSTRAINT_TOOLS}/bin:$PATH" 
PATH="${CONSTRAINT_TOOLS}/train-model:$PATH" 
PATH="${CONSTRAINT_TOOLS}/flask-app:${PATH}"

#######################################

## Define input variables
while getopts c:t: flag
do
	case "${flag}" in
		c) column=${OPTARG};;
		t) threshold=${OPTARG};;
	esac
done

info "Specifying depth cutoff to use for coverage filtering..."

## Determine column used for coverage filters
if [ ${column} == 5 ]
then
	column_name="over_1"
fi

if [ ${column} == 6 ]
then
	column_name="over_5"
fi

if [ ${column} == 7 ]
then
        column_name="over_10"
fi

if [ ${column} == 8 ]
then
        column_name="over_15"
fi

if [ ${column} == 9 ]
then
        column_name="over_20"
fi

if [ ${column} == 10 ]
then
        column_name="over_25"
fi

if [ ${column} == 11 ]
then
        column_name="over_30"
fi

if [ ${column} == 12 ]
then
        column_name="over_50"
fi

if [ ${column} == 13 ]
then
        column_name="over_100"
fi

if [ ${column} <= 4 | ${column} >= 14 ]
then
	echo "INVALID COLUMN NUMBER... PLEASE SPECIFY COLUMN VALUE BETWEEN 6 AND 13..."
fi

#######################################

## State coverage filter parameters
info "Filtering gnomad v3 coverage file to select for positions in which (${threshold} * 100)% of samples have a mean depth of value: ${column_name}X"

#######################################

## Define files to download
gnomad_url="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.1/vcf/genomes/gnomad.genomes.v3.1.1.sites.chr22.vcf.bgz"
gnomad_tbi_url="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.1/vcf/genomes/gnomad.genomes.v3.1.1.sites.chr22.vcf.bgz.tbi"
gnomad_coverage_v3="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.0.1/coverage/genomes/gnomad.genomes.r3.0.1.coverage.summary.tsv.bgz"

## Define directory to download files into 
gnomad_path="${CONSTRAINT_TOOLS}/data/gnomad/v3"

## Define inclusion/exclusion variant types
gnomad_inex_path="${CONSTRAINT_TOOLS}/data/gnomad/inclusion_exclusion_variants"
gnomad_inclusion=${gnomad_inex_path}/include_variants.txt
gnomad_exclusion=${gnomad_inex_path}/exclude_variants.txt

## Define reference genome
reference_genome="${CONSTRAINT_TOOLS}/data/reference/grch38/hg38.analysisSet.fa.gz"

## Define python script location to process gnomad files post-download
gnomad_scripts=${CONSTRAINT_TOOLS}/download-data/gnomad

## Define chr sizes input
chr_sizes="${CONSTRAINT_TOOLS}/data/chromosome-sizes/hg38.chrom.sizes.sorted"

## Define gnomad interval file
gnomad_intervals="${gnomad_path}/intervals/hg38.chrom.intervals.noheader.chr22"

mkdir --parents ${gnomad_path}

#######################################

## Define output names 
gnomad_variant="gnomad_v3_chr22.vcf.bgz"
gnomad_tbi="gnomad_v3_chr22.vcf.bgz.tbi"
gnomad_coverage="gnomad_v3_coverage.summary.tsv.bgz"
gnomad_coverage_filtered="gnomad_v3_coverage_${column_name}X_${threshold}.bed"
gnomad_vep_annotations="vep_annotations.gnomad_v3.txt"
gnomad_maf="gnomad_v3_chr22.maf"
chr_sizes_output="${CONSTRAINT_TOOLS}/data/gnomad/v3/intervals/hg38.chrom.intervals"

#######################################

info "Downloading gnomad's variant VCF file..."
#wget ${gnomad_url} --output-document=${gnomad_path}/${gnomad_variant}

info "Downloading gnomad's tbi file..."
#wget ${gnomad_tbi_url} --output-document=${gnomad_path}/${gnomad_tbi}

info "Downloading gnomad v3 coverage file..."
#wget ${gnomad_coverage_v3} --output-document=${gnomad_path}/${gnomad_coverage}

#######################################

info "Filtering gnomad v3 coverage file..."
#zcat ${gnomad_path}/${gnomad_coverage} --force | tail -n+2 | awk -v c=${column} -v t=${threshold} '$c>t {print $1}' | sed 's/:/\t/g' | awk '{print $1"\t"($2-1)"\t"$2}' | bedtools merge > ${gnomad_path}/${gnomad_coverage_filtered}

info "Identifying vep annotations for gnomad v3..."
#bcftools view -h ${gnomad_path}/${gnomad_variant} | grep "ID=vep" | grep "ID=vep" | tr ": " "\n" | grep Allele | sed 's/">//' | tr "|" "\n" > ${gnomad_path}/${gnomad_vep_annotations}

info "Segmenting chromosome sizes (hg38) for processing of gnomad v3 variant file..."
#python ${gnomad_scripts}/get_gnomad_v3_intervals.py --chr-sizes-file ${chr_sizes} --bin-num 100000 --output ${chr_sizes_output}  
#cat ${chr_sizes_output} | tail -n +2 > ${chr_sizes_output}.noheader

info "Processing gnomad v3 variant file..."
python ${gnomad_scripts}/process_gnomad_v3_variants.py --intervals ${gnomad_intervals} --gnomad_variant_file ${gnomad_path}/${gnomad_variant} --vep_annotation_file ${gnomad_path}/${gnomad_vep_annotations} --var_path ${gnomad_path}

info "Converting VCF variants to appropriate MAF format... Sorting and gzipping..."
#python ${gnomad_scripts}/vcf_to_maf.py --vcf ${gnomad_path}/gnomad_v3_variants.json --genome ${reference_genome} --kmer-size 3 --output ${gnomad_path}/${gnomad_maf}

info "Indexing gnomad v3 mutations..."
#cat ${gnomad_path}/${gnomad_maf} | bgzip > ${gnomad_path}/${gnomad_maf}.gz
#tabix \
	--skip-lines 1 \
	--sequence 1 \
	--begin 2 \
	--end 3 \
	--force \
	${gnomad_path}/${gnomad_maf}.gz
