#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH -o logs/constraint_tools_train.out

set -o errexit
set -o pipefail
set -o noclobber
set -o nounset

# set -o xtrace

## Define constraint tools working directory
CONSTRAINT_TOOLS=$1

source set-environment-variables.sh 

#######################################

# https://www.digitalocean.com/community/tutorials/how-to-read-and-set-environmental-and-shell-variables-on-linux

# https://stackoverflow.com/a/43476575/6674256
# need to export PYTHONPATH since it is not already in the environment: 
# `printenv | grep -w PYTHONPATH` returns zero output
export PYTHONPATH="${CONSTRAINT_TOOLS}/utilities:${CONSTRAINT_TOOLS}/predict-constraint"

# no need to export PATH since it is already in the environment: 
# `printenv | grep -w PATH` returns non-zero output 
PATH="${CONSTRAINT_TOOLS}/bin:$PATH"
PATH="${CONSTRAINT_TOOLS}/train-model:$PATH"

#######################################

info "Defining neutral region to train on..."

## Define neutral intervals to train on
neutral_regions="${CONSTRAINT_TOOLS}/dist/neutral-regions.bed.gz" 

## Remove regions less than 5 bp in length
#echo "" > ${CONSTRAINT_TOOLS}/dist/neutral-regions.filtered.bed
#zcat ${neutral_regions} | awk '{print $0"\t"$3-$2}' | awk '{ if ($4 >= 5) { print $1"\t"$2"\t"$3 }}' > ${CONSTRAINT_TOOLS}/dist/neutral-regions.filtered.bed

## Change neutral regions variable name
neutral_regions="/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/neutral-regions/neutral-regions.gnomadv3.filtered.bed"

## Define kmer size
kmer_size="5"

##################################################

info "Checking to see if potentially outdated model files exist in the output directory..."

## Define the model path for mutation probabilities output
model=${CONSTRAINT_TOOLS}/tests/train
mkdir --parents ${model}

## Check if the regional model directory is empty or not
#if [ $(find ${model} -empty -type f | wc -l) -eq 0 ]; then
	#info "regional model directory is empty... proceeding..."
#else 
	#info "regional model direcotry is not empty... removing potentially outdated contents..."
	#rm -rf ${model} 
	#mkdir --parents ${model}
#fi

##################################################

## Define number of jobs to run
job_num="500"

info "running constraint-tools train with ${job_num} jobs..."

## Define number of subjobs to run within each job --> a single subjob will train on a single neutral region
num_regions=($(cat ${neutral_regions} | wc -l))
subjob_num=$(printf "%.0f\n" $(echo ${num_regions} / ${job_num} + 0.5| bc -l))

##################################################

##################################################
##### Determine mutation pattern of interest #####
##################################################
cell_type="germline"

info "Looking at ${cell_type} constraint..."

## Get tumor sample information if the cell type of interest is "somatic"
if [ ${cell_type} == "somatic" ]; then

  # create and store a list of unique tumor barcodes, if such does not exist
  column_heading="Tumor_Sample_Barcode"
  if [[ ! -f ${mutations%.maf.gz}.${column_heading}.txt ]]; then
    info "Creating a list of unique values for the maf column: ${column_heading}"
    column_heading_index=$(fetch_column_heading_index ${mutations} ${column_heading})
    set +o errexit
    zcat ${mutations} |
      tail -n +2 | # lob off column headings
      cut -f ${column_heading_index} | # pull out column of interest
      # head -10000 | # debug
      sort | # required for uniq to work as expected
      uniq \
      > ${mutations%.maf.gz}.${column_heading}.txt
    set -o errexit
  fi

  number_samples=$(cat ${mutations%.maf.gz}.${column_heading}.txt | wc -l)
  info "Number tumors: ${number_tumors}\n"
fi

## Get the gnomad/topmed sample iformation if the cell type of interest is "germline"
if [ ${cell_type} == "germline" ]; then
  number_samples="76156" ## gnomad v3 WGS
  info "Number germline samples: ${number_samples}"
fi

##################################################

info "Running the job array... Training on ${subjob_num} neutral regions per job..."

## Define script to execute train
train_script="${CONSTRAINT_TOOLS}/train-model/train_on_region.sh"

## Generate job array sbatch script
#sbatch -W --array [1-${job_num}]%250 ${train_script} --CONSTRAINT-TOOLS-DIR ${CONSTRAINT_TOOLS} --kmer-size ${kmer_size} --cell-type ${cell_type} --model ${model} --neutral-regions-file ${neutral_regions} --subjob-num ${subjob_num} --number-samples ${number_samples}

##################################################

info "TODO (make this faster): Checking to see if all neutral regions could be trained on..."

list="ls ${model}"
#cat ${neutral_regions} | awk -v list="ls ${model}" '{print list " | grep " $1":"$2"-"$3}' | bash

##################################################

info "Concatenating model.json files obtained from each neutral region + calculating mutation probabilities..."

## Concatenate model.json files and calculate mutation probabilities

${CONSTRAINT_TOOLS}/train-model/calculate_mutation_probabilities \
  --final-model ${CONSTRAINT_TOOLS}/tests \
  --model-json-dir ${model}
