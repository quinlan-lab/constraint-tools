set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset 

source set-environment-variables.sh 

PATH="${CONSTRAINT_TOOLS}/experiments/germline-model/chen-et-al-2022:$PATH" 

# https://stackoverflow.com/a/43476575/6674256
# need to export PYTHONPATH since it is not already in the environment: 
# `printenv | grep -w PYTHONPATH` returns zero output
export PYTHONPATH="${CONSTRAINT_TOOLS}/utilities"

chen_mchale_windows="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/chen-mchale.bed"
snps="${CONSTRAINT_TOOLS_DATA}/gwas-catalog/gwas-catalog.tsv"
chen_mchale_windows_with_snp_counts="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/chen-mchale.snps.bed"

header-line () {
  echo -e "chromosome\tchen_start\tchen_end\tchen_zscore\tmchale_start\tmchale_end\tmchale_position\tmchale_N_bar\tmchale_N_observed\tmchale_K_bar\tmchale_K_observed\tmchale_M\toverlap_bps\tchen_snp_count"
}

get-chen-windows () { 
  cat ${chen_mchale_windows} | tail -n +2
}

get-unique-snp-positions () {
  get-snp-positions.py ${snps} \
    | sort --version-sort -k1,1 -k2,2 \
    | uniq 
}

add-snp-counts () {
  bedtools intersect \
      -a <(get-chen-windows) \
      -b <(get-unique-snp-positions) \
      -c
}

(
  header-line 
  add-snp-counts
) > ${chen_mchale_windows_with_snp_counts}  


