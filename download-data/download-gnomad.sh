#!/bin/sh
#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --nodes=2
#SBATCH --ntasks=16
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH -o download-gnomad.out

set -o errexit
set -o pipefail
# set -o noclobber
set -o xtrace
set -o nounset 

source download-data/set-environment-variables.sh 


## Define input variables
while getopts c:t: flag
do
	case "${flag}" in
		c) column=${OPTARG};;
		t) threshold=${OPTARG};;
	esac
done

## Determine column used for coverage filters
if [ $column == 5 ]
then
	column_name="over_1"
fi

if [ $column == 6 ]
then
	column_name="over_5"
fi

if [ $column == 7 ]
then
        column_name="over_10"
fi

if [ $column == 8 ]
then
        column_name="over_15"
fi

if [ $column == 9 ]
then
        column_name="over_20"
fi

if [ $column == 10 ]
then
        column_name="over_25"
fi

if [ $column == 11 ]
then
        column_name="over_30"
fi

if [ $column == 12 ]
then
        column_name="over_50"
fi

if [ $column == 13 ]
then
        column_name="over_100"
fi

if [ $column <= 4 | $column >= 14 ]
then
	echo "INVALID COLUMN NUMBER... PLEASE SPECIFY COLUMN VALUE BETWEEN 6 AND 13..."
fi

## State coverage filter parameters
info "Filtering gnomad v3 coverage file to select for positions in which ${threshold}% of samples have a mean depth of value: ${column_name}X"

## Define files to download
gnomad_url="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.1/vcf/genomes/gnomad.genomes.v3.1.1.sites.chr22.vcf.bgz"
gnomad_tbi_url="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.1/vcf/genomes/gnomad.genomes.v3.1.1.sites.chr22.vcf.bgz.tbi"
gnomad_coverage_v3="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.0.1/coverage/genomes/gnomad.genomes.r3.0.1.coverage.summary.tsv.bgz"

## Define directory to download files into 
gnomad_path="${CONSTRAINT_TOOLS}/data/gnomad/v3"

mkdir --parents ${gnomad_path}

## Define output names 
gnomad_variant="gnomad_v3_chr22.vcf.bgz"
gnomad_tbi="gnomad_v3_chr22.vcf.bgz.tbi"
gnomad_variant_filtered="gnomad_v3_chr22_filtered.vcf.bgz"
gnomad_coverage="gnomad_v3_coverage.summary.tsv.bgz"
gnomad_coverage_filtered="gnomad_v3_coverage_col${column}_${threshold}.bed"

info "Downloading gnomad's variant VCF file..."
wget ${gnomad_url} --output-document=${gnomad_path}/${gnomad_variant}

info "Downloading gnomad's tbi file..."
wget ${gnomad_tbi_url} --output-document=${gnomad_path}/${gnomad_tbi}

info "Downloading gnomad v3 coverage file..."
wget ${gnomad_coverage_v3} --output-document=${gnomad_path}/${gnomad_coverage}

info "Filtering gnomad v3 coverage file..."
zcat ${gnomad_path}/${gnomad_coverage} --force | tail -n+2 | awk -v c=${column} -v t=${threshold} '$c>t {print $1}' | sed 's/:/\t/g' | awk '{print $1"\t"($2-1)"\t"$2}' | bedtools merge > ${gnomad_path}/${gnomad_coverage_filtered}




