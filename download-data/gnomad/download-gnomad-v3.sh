#!/bin/sh
#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --nodes=4
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
      

#######################################

## Define files to download
gnomad_url="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.1/vcf/genomes/gnomad.genomes.v3.1.1.sites.chr22.vcf.bgz"
gnomad_tbi_url="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.1/vcf/genomes/gnomad.genomes.v3.1.1.sites.chr22.vcf.bgz.tbi"

## Define directory to download files into 
var_path="${CONSTRAINT_TOOLS}/data/gnomad/v3"

## Define reference genome
reference_genome="${CONSTRAINT_TOOLS}/data/reference/grch38/hg38.analysisSet.fa.gz"

## Define python script location to process gnomad files post-download
gnomad_scripts="${CONSTRAINT_TOOLS}/download-data/gnomad"

## Define chr sizes input
chr_sizes="${CONSTRAINT_TOOLS}/data/chromosome-sizes/hg38.chrom.sizes.sorted"

mkdir --parents ${var_path}

#######################################

## Define output names 
gnomad_variant="gnomad_v3_chr22.vcf.bgz"
gnomad_vep_annotations="vep_annotations.gnomad_v3.txt"
gnomad_maf="gnomad_v3_chr22.maf"
chr_intervals="${CONSTRAINT_TOOLS}/data/gnomad/v3/intervals/hg38.chrom.intervals"

#######################################

info "Downloading gnomad's variant VCF file..."
#wget ${gnomad_url} --output-document=${var_path}/vcf/${gnomad_variant}

info "Downloading gnomad's tbi file..."
#wget ${gnomad_tbi_url} --output-document=${var_path}/vcf/${gnomad_variant}.tbi

#######################################

info "Identifying vep annotations for gnomad v3..."
#bcftools view -h ${var_path}/${gnomad_variant} | grep "ID=vep" | tr ": " "\n" | grep Allele | sed 's/">//' | tr "|" "\n" > ${var_path}/vcf/${gnomad_vep_annotations}

info "Segmenting chromosome sizes (hg38) for processing of gnomad v3 variant file..."
#python ${gnomad_scripts}/get_gnomad_v3_intervals.py --chr-sizes-file ${chr_sizes} --bin-num 100000 --output ${chr_intervals}  
#cat ${chr_intervals} | tail -n +2 > ${chr_intervals}.noheader

info "Processing gnomad v3 variant file..."
python ${gnomad_scripts}/process_gnomad_v3_variants.py --intervals ${chr_intervals}.noheader.chr22 --gnomad-variant-file ${var_path}/vcf/${gnomad_variant} --vep-annotation-file ${var_path}/vcf/${gnomad_vep_annotations} --coverage-file ${var_path}/coverage/gnomad_v3_coverage.filtered.bed --var-path ${var_path}
