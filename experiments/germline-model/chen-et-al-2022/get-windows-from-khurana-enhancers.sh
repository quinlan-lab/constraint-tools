# data courtesy: Duo Xu
# git show --name-only c349b0bf7fd1ebb4433e4fd514d99ddd21466a67
ENHANCERS_STEM="/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/khurana/all-enhancers-with-network-features.hg38.sorted" 

get-windows-head () {
  jq -r '. | @tsv' "${ENHANCERS_STEM}.json" \
    | awk --assign OFS=$'\t' '{ print "window_chrom", "window_start", "window_end", $0 }'
}

# create windows by subtracting 0.5*WINDOW_SIZE from enhancer_midpoint and adding 0.5*WINDOW_SIZE to enhancer_midpoint
get-windows-tail () {
  cat "${ENHANCERS_STEM}.bed" \
    | get-regular-chromosomes \
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