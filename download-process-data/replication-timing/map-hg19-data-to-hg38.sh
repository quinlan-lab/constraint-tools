# set -o errexit
# set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset

source set-environment-variables.sh 

replication_timing_data_hg19="${CONSTRAINT_TOOLS_DATA}/replication-timing/iPSC_individual_level_data.hg19.txt.gz"
replication_timing_data_hg19_processed="${CONSTRAINT_TOOLS_DATA}/replication-timing/iPSC_individual_level_data.hg19.bed"
replication_timing_data_hg38="${CONSTRAINT_TOOLS_DATA}/replication-timing/iPSC_individual_level_data.hg38.bed" 

process-hg19-data () { 
  less ${replication_timing_data_hg19} \
    | python ${CONSTRAINT_TOOLS}/download-process-data/replication-timing/process_hg19_data.py \
    > ${replication_timing_data_hg19_processed}
}

map-hg19-data-to-hg38 () {
  # https://gist.github.com/brentp/894555/f23d1d6e0c988d6711acf2fe1a5bb930c3a19604
  # also see: https://genome.ucsc.edu/FAQ/FAQdownloads.html#liftOver
  # to explain errors, use "-errorHelp"
  wget --quiet --no-check-certificate -O - https://gist.github.com/raw/894555/lift.sh \
    | sh -s ${replication_timing_data_hg19_processed} hg19 hg38 "-bedPlus=3 -tab"
}

add-header () { 
  python ${CONSTRAINT_TOOLS}/download-process-data/replication-timing/add_header.py \
    ${replication_timing_data_hg19_processed} \
    ${replication_timing_data_hg19_processed}.hg38 \
  > ${replication_timing_data_hg38}
  rm \
    ${replication_timing_data_hg19_processed} \
    ${replication_timing_data_hg19_processed}.hg38
}

process-hg19-data
map-hg19-data-to-hg38
add-header
info "Wrote replication timing data in hg38 coordinates to" ${replication_timing_data_hg38}  

