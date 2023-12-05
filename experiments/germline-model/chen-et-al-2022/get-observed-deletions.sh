TOPMED="/scratch/ucgd/lustre-work/quinlan/u0055382/SVAFotate/supporting_data/TOPMed.GRCh38.bed.gz"

get-deletions-head () {
  less ${TOPMED} \
    | head -1 
}

# both het and homalt deletions: 
# filter out suspiciously large deletions:
get-deletions-tail () { 
  less ${TOPMED} \
    | tail -n +2 \
    | awk '{print "chr"$0}' \
    | awk '$5 == "DEL"' \
    | awk -v threshold=${SUSPICIOUS_DELETION_SIZE_THRESHOLD} '$4 < threshold'
}
