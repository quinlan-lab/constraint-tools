# experiments/germline-model/chen-et-al-2022/SNV_plus_SV_model.4.ipynb
ENHANCERS="/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/chen-et-al-2023-published-version/41586_2023_6045_MOESM4_ESM/Supplementary_Data_6_ESM-intersect-all-observed-deletions-with-labels.bed" 

get-windows-head () { 
  set +o errexit
  cat ${ENHANCERS} \
    | head -1 \
    | awk --assign OFS=$'\t' '{ print "window_chrom", "window_start", "window_end", $0 }'
  set -o errexit
}

# TODO: 
# generalize this code to use "bedtools flank"
get-windows-tail () {
  cat ${ENHANCERS} \
    | tail -n +2 \
    | awk  \
      --assign OFS=$'\t' \
      --assign window_size=${WINDOW_SIZE} \
      '{ 
        midpoint = 0.5*($2+$3)
        start = (midpoint-0.5*window_size > 0 ? midpoint-0.5*window_size : 0)
        end = midpoint+0.5*window_size
        printf "%s\t%d\t%d\t%s\n", $1, start, end, $0
      }'
}

