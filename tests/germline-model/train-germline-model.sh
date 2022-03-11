set -o errexit
set -o pipefail
set -o nounset

number_of_jobs="500"

CONSTRAINT_TOOLS=$PWD
PATH=${CONSTRAINT_TOOLS}/utilities:$PATH 

training_durations_filename="tests/germline-model/train-durations.${number_of_jobs}-jobs.tsv"
echo -e "number_of_neutral_regions\ttraining_duration" > ${training_durations_filename} 

for number_of_neutral_regions in "10000" "100000" "1000000"; do
  info "Number of (unfiltered) neutral regions:" ${number_of_neutral_regions}
  time_output_file="tests/germline-model/${number_of_neutral_regions}.time"
  /usr/bin/time --verbose --output ${time_output_file} bash train-germline-model.sh \
    --constraint-tools-directory ${CONSTRAINT_TOOLS} \
    --number-of-jobs ${number_of_jobs} \
    --number-of-neutral-regions ${number_of_neutral_regions} \
    --work-directory "work-train-germline-model-${number_of_neutral_regions}" \
    --work-directory-should-be-clean "true"
  info "time output file:" ${time_output_file}
  training_duration=$(grep "Elapsed (wall clock) time" ${time_output_file} | awk '{ print $NF }')
  info "Training duration:" ${training_duration}
  echo -e "${number_of_neutral_regions}\t${training_duration} (h:mm:ss or m:ss)" >> ${training_durations_filename}
  echo -e "\n***************************************************"
done 

