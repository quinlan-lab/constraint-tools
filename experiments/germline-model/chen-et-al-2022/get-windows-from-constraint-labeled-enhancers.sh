# experiments/germline-model/chen-et-al-2022/SNV_plus_SV_model.4.ipynb
ENHANCERS="/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/chen-et-al-2023-published-version/41586_2023_6045_MOESM4_ESM/Supplementary_Data_6_ESM-intersect-all-observed-deletions-with-labels.bed" 

get-windows-head () { 
  set +o errexit
  cat ${ENHANCERS} \
    | head -1 \
    | awk --assign OFS=$'\t' '{ print "window_chrom", "window_start", "window_end", $0 }'
  set -o errexit
}

# create windows by subtracting 0.5*WINDOW_SIZE from enhancer_midpoint and adding 0.5*WINDOW_SIZE to enhancer_midpoint
get-windows-tail () {
  cat ${ENHANCERS} \
    | tail -n +2 \
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

