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

## State coverage filter parameters
info "Filtering gnomad v3 coverage file to select for positions in which ${threshold}% of samples have a mean depth greater than the value specified by column number ${column}"

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
gnomad_coverage="gnomad_v3_coverage.summary.tsv.bgz"
gnomad_coverage_filtered="gnomad_v3_coverage_${threshold}.bed"

info "Downloading gnomad's variant VCF file..."
#wget ${gnomad_url} --output-document=${gnomad_path}/${gnomad_variant}

info "Downloading gnomad's tbi file..."
#wget ${gnomad_tbi_url} --output-document=${gnomad_path}/${gnomad_tbi}

info "Downloading gnomad v3 coverage file..."
#wget ${gnomad_coverage_v3} --output-document=${gnomad_path}/${gnomad_coverage}

info "Filtering gnomad v3 coverage file..."
zcat ${gnomad_path}/${gnomad_coverage} --force | tail -n+2 | awk -v c=${column} -v t=${threshold} '$c>t {print $1}' | sed 's/:/\t/g' | awk '{print $1"\t"($2-1)"\t"$2}' | bedtools merge > ${gnomad_path}/${gnomad_coverage_filtered}


