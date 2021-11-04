#!/bin/sh
#SBATCH --time=20:00:00
#SBATCH --nodes=4
#SBATCH --ntasks=16
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH -o logs/predict_constraint/predict-constraint-main.out

set -o errexit
set -o pipefail
# set -o noclobber
set -o xtrace
set -o nounset

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

## Define main predict constraint directory (data)
prediction_constraint_dir="${CONSTRAINT_TOOLS}/predict-constraint"

## Get file with transcripts to predict constraint on
transcripts="${prediction_constraint_dir}/prediction_intervals/transcripts.txt"

info "Defining job array parameters..."

## Define number of transcripts to predict constraint on
transcript_num=($(cat ${transcripts}| wc -l))

## Define size of job array
job_array_size="500"

## Define parameters for job array
subjob_num=$(echo "(${transcript_num}+1)/500"+1 | bc)

## Submit the job array to predict constraint
log_dir="${CONSTRAINT_TOOLS}/logs/predict_constraint/prediction_logs/"
mkdir --parents ${log_dir}
log="${log_dir}/predict_constraint-%a.out"

info "Predicting constraint for each transcript in ${transcripts}..."

## Define main directory for predict constraint scripts
predict_constraint_scripts="${CONSTRAINT_TOOLS}/predict-constraint"

## Define json model
model_filename="${CONSTRAINT_TOOLS}/tests/model.json"

## Define the coverage file
coverage_filename="${CONSTRAINT_TOOLS}/dist/covered_sites_wgs.sorted.bed"

## Define output directory for predict constraint
output_dir="${CONSTRAINT_TOOLS}/predict-constraint/intermediate_files/transcript_predict_constraint"
mkdir --parents ${output_dir}

sbatch -W --array [1-${job_array_size}]%50 --output ${log} ${predict_constraint_scripts}/predict_constraint.sh --transcripts ${transcripts} --model-filename ${model_filename} --coverage-filename ${coverage_filename} --output-path ${output_dir}

info "Concatenating all of the prediction files..."

## Concatenate all of the prediction files
predict_constraint_file="${CONSTRAINT_TOOLS}/predict-constraint/final/predict_constraint_final.bed"
echo -n "" > ${predict_constraint_file}
cat ${output_dir}/* >> ${predict_constraint_file}

