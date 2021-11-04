#!/bin/sh
#SBATCH --time=3:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH -o logs/predict_constraint/get_prediction_intervals.out

set -o errexit
set -o pipefail
# set -o noclobber
set -o xtrace
set -o nounset

source download-data/set-environment-variables.sh 

## Define gtf file
gtf_file="${CONSTRAINT_TOOLS}/data/genes/Homo_sapiens.GRCh38.104.sorted.gtf.gz"

## Define interval files and directory
prediction_intervals_dir="${CONSTRAINT_TOOLS}/predict-constraint/prediction_intervals"
predict_constraint_intervals="${prediction_intervals_dir}/predict_constraint_intervals"

## Filter gtf file for protein coding exons
echo -n "" > ${predict_constraint_intervals}.bed
zcat ${gtf_file} | grep -w "gene_biotype \"protein_coding\"" | grep -w "CDS" | cut -f 1,4,5,9 > ${predict_constraint_intervals}.bed 

## Get list of transcripts
transcripts="${prediction_intervals_dir}/transcripts.txt"
cat ${predict_constraint_intervals}.bed | cut -f 4 | sed 's/; /\t/g' | cut -f 3 | sed 's/"//g' | sed 's/[^ ]* //' | sort | uniq > ${transcripts}

##### TODO: Make the above command more robust: 
## Ensure that the appropriate column/value is being extracted from the gtf file
