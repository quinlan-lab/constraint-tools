CHROMOSOME_SIZES="${CONSTRAINT_TOOLS_DATA}/reference/grch38/chromosome-sizes/hg38.chrom.sizes.sorted"

labeled_enhancers_filename="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/labeled-enhancers.gnocchi.GC.bed"

get-labeled-enhancers-tail () { 
  tail -n +2 "${labeled_enhancers_filename}" \
    | cut -f1-3
}

get-GC-windows-head () {
  echo -e "GC_window_chrom\tGC_window_start\tGC_window_end\tlabeled_enhancer_chrom\tlabeled_enhancer_start\tlabeled_enhancer_end"  
}

# create GC-windows by subtracting 0.5*WINDOW_SIZE from labeled_enhancer_midpoint and adding 0.5*WINDOW_SIZE to labeled_enhancer_midpoint
# note that "set -o nounset" does not catch unbound variables in "$((WINDOW_SIZE/2))"
get-GC-windows-tail () {
  get-labeled-enhancers-tail \
    | awk  \
      --assign OFS=$'\t' \
      '{ 
        midpoint = 0.5*($2+$3)
        start = midpoint-1 
        end = midpoint
        printf "%s\t%d\t%d\t%s\n", $1, start, end, $0
      }' \
    | bedtools slop -i - -g ${CHROMOSOME_SIZES} -b $((WINDOW_SIZE/2))
}
