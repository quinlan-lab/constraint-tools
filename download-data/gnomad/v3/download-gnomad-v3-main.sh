#!/bin/sh
#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH -o logs/download-gnomad-v3-main.out


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

## Define location of gnomad v3 download scripts
gnomad_scripts="${CONSTRAINT_TOOLS}/download-data/gnomad/v3"

#######################################                                                                  
var_path="/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/gnomad/v3/vcf"
mkdir --parents ${var_path}

#for chromosome in $(seq 1 22) X Y
for chromosome in Y
do 
	chromosome="chr${chromosome}"
	gnomad_variant="gnomad_v3_${chromosome}"

	info "Downloading and extracting VCF info fields from gnomad v3's ${chromosome} vcf file..."

	info "Making directory for ${chromosome} log files..."
	log_dir="${CONSTRAINT_TOOLS}/logs/gnomad_v3_variants/${chromosome}"
	mkdir --parents ${log_dir}
	
	info "Executing download and processing script..."
	sbatch -o ${log_dir}/${chromosome}_download_gnomad_v3_vcf.out ${gnomad_scripts}/download-gnomad-v3-vcf.sh --chromosome ${chromosome} --var-path ${var_path} --gnomad-variant ${gnomad_variant}
done

