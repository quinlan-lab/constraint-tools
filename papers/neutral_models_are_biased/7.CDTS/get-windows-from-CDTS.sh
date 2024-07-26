CHROMOSOME_SIZES="${CONSTRAINT_TOOLS_DATA}/reference/grch38/chromosome-sizes/hg38.chrom.sizes.sorted"

CDTS_filename="${CONSTRAINT_TOOLS_DATA}/CDTS/CDTS.gnomAD.hg38.noncoding.enhancer.bed"

get-CDTS-windows-tail () { 
  tail -n +2 "${CDTS_filename}" \
    | cut -f1-3
}

get-GC-windows-head () {
  echo -e "GC_window_chrom\tGC_window_start\tGC_window_end\tCDTS_window_chrom\tCDTS_window_start\tCDTS_window_end"  
}

# create GC-windows by subtracting 0.5*WINDOW_SIZE from CDTS_window_midpoint and adding 0.5*WINDOW_SIZE to CDTS_window_midpoint
# note that "set -o nounset" does not catch unbound variables in "$((WINDOW_SIZE/2))"
get-GC-windows-tail () {
  get-CDTS-windows-tail \
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