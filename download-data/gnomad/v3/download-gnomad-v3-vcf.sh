#!/bin/sh
#!/bin/bash
#SBATCH --time=71:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw

## Define input variables
while [[ "$1" =~ ^- ]]; do
  case $1 in
    --chromosome ) shift; [[ ! $1 =~ ^- ]] && chromosome=$1;;
    --var-path ) shift; [[ ! $1 =~ ^- ]] && var_path=$1;;
    --gnomad-variant ) shift; [[ ! $1 =~ ^- ]] && gnomad_variant=$1;;
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

# https://stackoverflow.com/a/43476575/6674256
# need to export PYTHONPATH since it is not already in the environment:
# `printenv | grep -w PYTHONPATH` returns zero output
export PYTHONPATH="${CONSTRAINT_TOOLS}/utilities:${CONSTRAINT_TOOLS}/predict-constraint"

# no need to export PATH since it is already in the environment:
# `printenv | grep -w PATH` returns non-zero output
PATH="${CONSTRAINT_TOOLS}/bin:$PATH"

#######################################

info "Downloading gnomad v3's ${chromosome} vcf file..."

## Define files to download
gnomad_url="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.1/vcf/genomes/gnomad.genomes.v3.1.1.sites.${chromosome}.vcf.bgz"
gnomad_tbi_url="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.1/vcf/genomes/gnomad.genomes.v3.1.1.sites.${chromosome}.vcf.bgz.tbi"

info "Downloading gnomad's variant VCF file..."
wget ${gnomad_url} --output-document=${var_path}/${gnomad_variant}.vcf.gz

info "Downloading gnomad's tbi file..."
wget ${gnomad_tbi_url} --output-document=${var_path}/${gnomad_variant}vcf.gz.tbi

info "Done..."
