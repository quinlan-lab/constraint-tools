#!/bin/sh
#SBATCH --time=23:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw

## Define input variables
while [[ "$1" =~ ^- ]]; do
  case $1 in
    --transcripts ) shift; [[ ! $1 =~ ^- ]] && transcripts=$1;;
    --model-filename ) shift; [[ ! $1 =~ ^- ]] && model_filename=$1;;
    --coverage-filename ) shift; [[ ! $1 =~ ^- ]] && coverage_filename=$1;;
    --output-path ) shift; [[ ! $1 =~ ^- ]] && output_path=$1;;
    *) error "$0: $1 is an invalid flag"; exit 1;;
  esac
  shift
done


set -o errexit
set -o pipefail
# set -o noclobber
set -o xtrace
set -o nounset

source download-data/set-environment-variables.sh

# no need to export PATH since it is already in the environment:
# `printenv | grep -w PATH` returns non-zero output
PATH="${CONSTRAINT_TOOLS}/bin:$PATH"

####################

info "Setting up subjobs for job array number: ${SLURM_ARRAY_TASK_ID}..."

## The job array size will be the number of lines in the dile divided by the number of lines chosen below
start=${SLURM_ARRAY_TASK_ID}
numlines=5
stop=$((${SLURM_ARRAY_TASK_ID}*numlines))
start="$((${stop} - $((numlines -1))))"

echo "start=${start}"
echo "stop=${stop}"

for (( line = ${start}; line <= ${stop}; line++ ))
do

        ########################################

	info "Defining transcript of interest..."

	## Get the transcript of interest
	transcript=$(sed -n ${line}p ${transcripts})

        ########################################

	info "Identifying protein coding exons for this transcript..."
	## Define directory for intermediate files
	transcript_intermediate_files_dir="${CONSTRAINT_TOOLS}/predict-constraint/intermediate_files/intervals/${transcript}"
	mkdir --parents ${transcript_intermediate_files_dir} 

	## Define transcript-specific intermediate file name
	transcript_bed="${transcript_intermediate_files_dir}/${transcript}.sorted.bed"

	## Subset genome-wide cds exons for transcript of interest
	all_exons="${CONSTRAINT_TOOLS}/predict-constraint/prediction_intervals/predict_constraint_intervals.bed"
	cat ${all_exons} | grep -w "transcript_id \"${transcript}\"" | awk -v OFS='\t' '{print "chr"$1,$2,$3}' | sort -k1,1 -k2,2n > ${transcript_bed}	

	########################################

	info "Identifying portions of ${transcript}'s exons that are sufficiently covered..."

	transcript_covered_bed="${transcript_intermediate_files_dir}/${transcript}.covered.bed"
	cat ${transcript_bed} | bedtools intersect -a - -b ${coverage_filename} -sorted > ${transcript_covered_bed}

        ########################################

        info "Identifying portions of ${transcript}'s exons that are INsufficiently covered..."

	transcript_notcovered_bed="${transcript_intermediate_files_dir}/${transcript}.notcovered.bed"
	bedtools subtract -a ${transcript_bed} -b ${transcript_covered_bed} > ${transcript_notcovered_bed}

        ########################################

	info "Combining ${transcript}'s covered and not covered regions..."
	
	## Add covered and notcovered tages
	cat ${transcript_covered_bed} | awk '{print $0, "\t", "covered"}' > ${transcript_covered_bed}.tagged
	cat ${transcript_notcovered_bed} | awk '{print $0, "\t", "not_covered"}' > ${transcript_notcovered_bed}.tagged

	transcript_bed_final="${transcript_intermediate_files_dir}/${transcript}.sorted.final.bed"
	cat ${transcript_covered_bed}.tagged ${transcript_notcovered_bed}.tagged | sort -k1,1 -k2,2n > ${transcript_bed_final}
	
        ########################################

	info "Predicting constraint on ${transcript}..."

	## Run predict-constraint python script 
	predict_constraint_scripts="${CONSTRAINT_TOOLS}/predict-constraint"
	window_size="101" ## Has to be odd
	window_stride="50"
	python ${predict_constraint_scripts}/predict_constraint.py --transcript ${transcript} --transcript-intervals ${transcript_bed_final} --window-size ${window_size} --window-stride ${window_stride} --model-filename ${model_filename} --output-path ${output_path}

	#rm ${transcript_bed}
done




