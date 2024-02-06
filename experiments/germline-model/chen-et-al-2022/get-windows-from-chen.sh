CHROMOSOME_SIZES="${CONSTRAINT_TOOLS_DATA}/reference/grch38/chromosome-sizes/hg38.chrom.sizes.sorted"

get-chen-windows () { 
  # chen scores for 1kb windows: 
  local filename="${CHEN_DATA_DIRECTORY}/${CHEN_FILE_STEM}.bed"
  cat "${filename}"
}

get-windows-head () { 
  echo -e "window_chrom\twindow_start\twindow_end\tchen_chrom\tchen_start\tchen_end\tchen_score"
}

# create windows by subtracting 0.5*WINDOW_SIZE from chen_window_midpoint and adding 0.5*WINDOW_SIZE to chen_window_midpoint
# note that "set -o nounset" does not catch unbound variables in "$((WINDOW_SIZE/2))"
get-windows-tail () {
  get-chen-windows \
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

