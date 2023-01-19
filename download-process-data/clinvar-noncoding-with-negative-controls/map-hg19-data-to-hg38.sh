set -o errexit
set -o pipefail
# set -o xtrace
set -o nounset

source download-process-data/set-environment-variables.sh 

data_hg19="${CONSTRAINT_TOOLS_DATA}/clinvar-noncoding-with-negative-controls/Clinvar_nc_snv_pathogenic_177Pos_177Neg.hg19.vcf"
liftover_unstable_positions="${CONSTRAINT_TOOLS_DATA}/clinvar-noncoding-with-negative-controls/FASTA_BED.ALL_GRCh37.novel_CUPs.bed"
data_hg19_processed="${CONSTRAINT_TOOLS_DATA}/clinvar-noncoding-with-negative-controls/Clinvar_nc_snv_pathogenic.hg19.bed"
data_hg38="${CONSTRAINT_TOOLS_DATA}/clinvar-noncoding-with-negative-controls/Clinvar_nc_snv_pathogenic.hg38.bed" 

process-hg19-data () { 
  less ${data_hg19} \
    | python ${CONSTRAINT_TOOLS}/download-process-data/clinvar-noncoding-with-negative-controls/process_hg19_data.py \
    | bedtools intersect -v -a - -b ${liftover_unstable_positions} \
    > ${data_hg19_processed}
}

map-hg19-data-to-hg38 () {
  # https://gist.github.com/brentp/894555/f23d1d6e0c988d6711acf2fe1a5bb930c3a19604
  # also see: https://genome.ucsc.edu/FAQ/FAQdownloads.html#liftOver
  # to explain errors, use "-errorHelp"
  wget --quiet --no-check-certificate -O - https://gist.github.com/raw/894555/lift.sh \
    | sh -s ${data_hg19_processed} hg19 hg38 "-bedPlus=3 -tab"
}

add-header () { 
  python ${CONSTRAINT_TOOLS}/download-process-data/clinvar-noncoding-with-negative-controls/add_header.py \
    ${data_hg19_processed}.hg38 \
  > ${data_hg38}
  rm \
    ${data_hg19_processed} \
    ${data_hg19_processed}.hg38
}

process-hg19-data
map-hg19-data-to-hg38
add-header
info "Wrote data in hg38 coordinates to" ${data_hg38}  

