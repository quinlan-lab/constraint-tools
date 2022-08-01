set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset 

source set-environment-variables.sh 

chen_windows="${CONSTRAINT_TOOLS_DATA}/chen-et-al-2022/Supplementary_Datasets/Supplementary_Data_2.bed"
mchale_windows="${CONSTRAINT_TOOLS_DATA}/genome-wide-predictions/predict-germline-grch38.windowSize-1001.bed.gz"
chen_mchale_windows="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/chen-mchale.bed"

header-line () {
  echo -e "chromosome\tchen_start\tchen_end\tchen_zscore\tmchale_start\tmchale_end\tmchale_position\tmchale_N_bar\tmchale_N_observed\tmchale_K_bar\tmchale_K_observed\tmchale_M\toverlap_bps"
}

overlapping-zscores () {
  local overlap_fraction="0.9" # should be >0.5 to ensure that each McHale window matches a single Chen window, and vice versa

  bedtools intersect \
      -a <(less ${chen_windows} | tail -n +2) \
      -b <(less ${mchale_windows} | tail -n +2) \
      -f ${overlap_fraction} \
      -r \
      -wo \
    | cut -f 1-4,6-14
}

(
  header-line 
  overlapping-zscores
) > ${chen_mchale_windows}  


