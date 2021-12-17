#!/bin/sh
#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw

## Define input variables
while [[ "$1" =~ ^- ]]; do
  case $1 in
    --chromosome ) shift; [[ ! $1 =~ ^- ]] && chromosome=$1;;
    --chr-interval-file ) shift; [[ ! $1 =~ ^- ]] && chr_interval_file=$1;;
    --gnomad-variant-file ) shift; [[ ! $1 =~ ^- ]] && gnomad_variant_file=$1;;
    --vep-annotation-file ) shift; [[ ! $1 =~ ^- ]] && vep_annotation_file=$1;;
    --var-path ) shift; [[ ! $1 =~ ^- ]] && var_path=$1;;
    *) error "$0: $1 is an invalid flag"; exit 1;;
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

## Make directory to store intermediate files
intermediate_files_dir="${var_path}/intermediate_files/${chromosome}"

## The job array size will be the number of lines in the dile divided by the number of lines chosen below
start=${SLURM_ARRAY_TASK_ID}
numlines=5
stop=$((${SLURM_ARRAY_TASK_ID}*numlines))
start="$((${stop} - $((numlines -1))))"

echo "start=${start}"
echo "stop=${stop}"

for (( line = ${start}; line <= ${stop}; line++ ))
do
	## Get the interval
	interval=$(sed -n ${line}p ${chr_interval_file} | awk '{print $1":"$2"-"$3}')
	
	info "Processing variants within ${interval}..."

	## Process variants within each interval	
	gnomad_scripts="${CONSTRAINT_TOOLS}/download-data/gnomad/v3"
	python ${gnomad_scripts}/process_gnomad_v3_variants.py --interval ${interval} --gnomad-variant-file ${gnomad_variant_file} --vep-annotation-file ${vep_annotation_file} --var-path ${var_path}
done
