#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

number_of_jobs=$1
number_of_trustworthy_noncoding_regions_list=$2
number_of_trustworthy_noncoding_regions_list=${number_of_trustworthy_noncoding_regions_list//,/ }

CONSTRAINT_TOOLS=$PWD
PATH=${CONSTRAINT_TOOLS}/utilities:$PATH 

training_durations_filename="tests/germline-model/train-durations.${number_of_jobs}-jobs.tsv"
echo -e "number_of_trustworthy_noncoding_regions\ttraining_duration" > "${training_durations_filename}" 

for number_of_trustworthy_noncoding_regions in $number_of_trustworthy_noncoding_regions_list; do
  info "Number of (unfiltered) trustworthy noncoding regions:" "${number_of_trustworthy_noncoding_regions}"
  time_output_file="tests/germline-model/${number_of_trustworthy_noncoding_regions}.time"
  /usr/bin/time --verbose --output "${time_output_file}" bash tests/germline-model/train-germline-model-core.sh \
    --constraint-tools-directory "${CONSTRAINT_TOOLS}" \
    --number-of-jobs "${number_of_jobs}" \
    --number-of-trustworthy-noncoding-regions "${number_of_trustworthy_noncoding_regions}" \
    --work-directory "work-train-germline-model-${number_of_trustworthy_noncoding_regions}" \
    --work-directory-should-be-clean "true"
  info "time output file:" "${time_output_file}"
  training_duration=$(grep "Elapsed (wall clock) time" "${time_output_file}" | awk '{ print $NF }')
  info "Training duration:" "${training_duration}"
  echo -e "${number_of_trustworthy_noncoding_regions}\t${training_duration} (h:mm:ss or m:ss)" >> "${training_durations_filename}"
  echo -e "\n***************************************************"
done 

