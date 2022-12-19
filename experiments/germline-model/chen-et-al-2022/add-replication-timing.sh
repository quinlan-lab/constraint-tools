set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset 

source set-environment-variables.sh 

chen_windows="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/chen-mchale-enhancer-exon.bed"
replication_timing_windows="${CONSTRAINT_TOOLS_DATA}/replication-timing/iPSC_individual_level_data.hg38.bed"
chen_windows_with_replication_timing_data="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/chen-mchale-enhancer-exon-replicationTiming.bed"

header-line () {
  echo -e "chromosome\tstart\tend\tposition\tN_bar\tN_observed\tK_bar\tK_observed\tM\tchen zscore\tenhancer overlap\tmerged_exon overlap\twindow overlaps enhancer\twindow overlaps merged_exon\twindow overlaps (enhancer, merged_exon)\tnegative chen zscore\tsingleton proportion\tsingleton proportion error (ind. sites)\tnegative N_bar\tchromosome_RT\tstart_RT\tend_RT\treplication timing per individual\tchen-replicationTiming window overlap"
}

add-replication-timing () {
  local min_overlap_as_a_fraction_of_chen_window="0.9" # should be >0.5 to ensure that each Chen window matches at most a single replication-timing window

  bedtools intersect \
      -a <(less ${chen_windows} | tail -n +2) \
      -b <(less ${replication_timing_windows} | tail -n +2) \
      -f ${min_overlap_as_a_fraction_of_chen_window} \
      -wo 
}

(
  header-line 
  add-replication-timing 
) > ${chen_windows_with_replication_timing_data}  

info "wrote to:" ${chen_windows_with_replication_timing_data}  


