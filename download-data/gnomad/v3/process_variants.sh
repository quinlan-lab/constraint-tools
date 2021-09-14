#!/bin/sh
#!/bin/bash

## Define input variables
while [[ "$1" =~ ^- ]]; do
  case $1 in
    --chromosome ) shift; [[ ! $1 =~ ^- ]] && chromosome=$1;;
    --interval ) shift; [[ ! $1 =~ ^- ]] && interval=$1;;
    --gnomad-variant-file ) shift; [[ ! $1 =~ ^- ]] && gnomad_variant_file=$1;;
    --vep-annotation-file ) shift; [[ ! $1 =~ ^- ]] && vep_annotation_file=$1;;
    --var-path ) shift; [[ ! $1 =~ ^- ]] && var_path=$1;;
    *) error "$0: $1 is an invalid flag"; exit 1;;
  esac
  shift
done

#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH -o logs/gnomad_v3_variants/${chromosome}/download-gnomad-v3-%j-${chromosome}.out

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
gnomad_scripts="${CONSTRAINT_TOOLS}/download-data/gnomad/v3"
python ${gnomad_scripts}/process_gnomad_v3_variants.py --interval ${interval} --gnomad-variant-file ${gnomad_variant_file} --vep-annotation-file ${vep_annotation_file} --var-path ${var_path}

