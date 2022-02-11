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

number_tumors=$(cat ${mutations%.maf.gz}.${column_heading}.txt | wc -l)
info "Number tumors: ${number_tumors}\n"



