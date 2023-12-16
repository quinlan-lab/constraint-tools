# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5467550/
# data courtesy: Tom Nicholas 
# https://docs.google.com/presentation/d/1o1rOcvnj69F9-qfXK0r19UKLupoivBvBlzBboM9u2dQ/edit#slide=id.g22aad12e963_0_24
ENHANCERS="/scratch/ucgd/lustre-work/quinlan/u0055382/genome_reference/genehancer/genehancer_GRCh38.bed.gz" 

get-windows-head () { 
  set +o errexit
  zcat ${ENHANCERS} \
    | head -1 \
    | awk --assign OFS=$'\t' '{ print "window_chrom", "window_start", "window_end", $0 }'
  set -o errexit
}

# TODO: 
# generalize this code to use "bedtools flank"
get-windows-tail () {
  zcat ${ENHANCERS} \
    | tail -n +2 \
    | awk  \
      --assign OFS=$'\t' \
      --assign window_size=${WINDOW_SIZE} \
      '{ 
        midpoint = 0.5*($2+$3)
        start = (midpoint-0.5*window_size > 0 ? midpoint-0.5*window_size : 0)
        end = midpoint+0.5*window_size
        printf "chr%s\t%d\t%d\tchr%s\n", $1, start, end, $0
      }'
}
