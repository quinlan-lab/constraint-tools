#!/bin/bash

#SBATCH --time=71:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH -o logs/train/train-model-%a.out

while [[ "$1" =~ ^- ]]; do 
  case $1 in
    --CONSTRAINT-TOOLS-DIR ) shift; [[ ! $1 =~ ^- ]] && CONSTRAINT_TOOLS=$1;;
    --kmer-size ) shift; [[ ! $1 =~ ^- ]] && kmer_size=$1;;
    --cell-type ) shift; [[ ! $1 =~ ^- ]] && cell_type=$1;;
    --model ) shift; [[ ! $1 =~ ^- ]] && model=$1;;
    --neutral-regions-file ) shift; [[ ! $1 =~ ^- ]] && neutral_regions_file=$1;;
    --subjob-num ) shift; [[ ! $1 =~ ^- ]] && subjob_num=$1;;
    --number-samples ) shift; [[ ! $1 =~ ^- ]] && number_samples=$1;;
    *) error "$0: $1 is an invalid flag"; exit 1;;
  esac 
  shift
done

set -o errexit
set -o pipefail
set -o nounset
# set -o noclobber
# set -o xtrace

#######################################

## Define mutation file and reference genome based on analysis type (i.e. somatic vs germline)
if [ ${cell_type} == "somatic" ]; then
	mutations="${CONSTRAINT_TOOLS}/data/icgc/mutations.sorted.maf.gz"
	genome="${CONSTRAINT_TOOLS}/data/reference/grch37/genome.fa.gz"

elif [ ${cell_type} == "germline" ]; then
	mutations="/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/gnomad/v3/gnomad_v3_variants.sorted.bed.gz"
	genome="${CONSTRAINT_TOOLS}/data/reference/grch38/hg38.analysisSet.fa.gz"

else 
	info "PLEASE SUPPLY \"germline\" OR \"somatic\" as input for the cell_type variable..."
fi

#######################################

## The job array size will be the number of lines in the file divided by the number of lines chosen below
start=${SLURM_ARRAY_TASK_ID}
numlines=${subjob_num}
stop=$((${SLURM_ARRAY_TASK_ID}*numlines))
start="$((${stop} - $((numlines -1))))"

echo "start=${start}"
echo "stop=${stop}"

#######################################

for (( line = ${start}; line <= ${stop}; line++))
do
	## Determine region to train on
	region=$(cat ${neutral_regions_file} | sed -n ${line}p | awk '{print $1":"$2"-"$3}')

	echo ${region}
	## Skip if there is no region, this might occur when looking at the last subjob of the 500th iteration of the job array
	if [ -z ${region} ]; then
		continue
	fi
	
	${CONSTRAINT_TOOLS}/constraint-tools train \
  	  --genome ${genome} \
  	  --mutations ${mutations} \
  	  --kmer-size ${kmer_size} \
  	  --region ${region} \
  	  --model ${model} \
 	  --number-samples ${number_samples}
done
