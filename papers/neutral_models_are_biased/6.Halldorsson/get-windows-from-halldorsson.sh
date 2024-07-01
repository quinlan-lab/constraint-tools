CHROMOSOME_SIZES="${CONSTRAINT_TOOLS_DATA}/reference/grch38/chromosome-sizes/hg38.chrom.sizes.sorted"

get-halldorsson-windows-tail () { 
  local filename="${CONSTRAINT_TOOLS_DATA}/depletion_rank_scores/41586_2022_4965_MOESM3_ESM"
  tail -n +2 "${filename}"
}

get-windows-head () { 
  echo -e "window_chrom\twindow_start\twindow_end\thalldorsson_chrom\thalldorsson_start\thalldorson_end\thalldorsson_score"
}

# create windows by subtracting 0.5*WINDOW_SIZE from halldorsson_window_midpoint and adding 0.5*WINDOW_SIZE to halldorsson_window_midpoint
# note that "set -o nounset" does not catch unbound variables in "$((WINDOW_SIZE/2))"
get-windows-tail () {
  get-halldorsson-windows-tail \
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

