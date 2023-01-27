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
                                              
chen_mchale_windows="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/chen-mchale-enhancer-exon.bed"
clinvar_snvs="${CONSTRAINT_TOOLS_DATA}/clinvar-noncoding-with-negative-controls/Clinvar_nc_snv_pathogenic.hg38.bed"
chen_mchale_windows_with_clinvar_snv_counts="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/chen-mchale-enhancer-exon-clinvar.bed"

header-line () {
  set +o errexit
  echo -e "$(head -1 ${chen_mchale_windows})\tClinVar SNV count"
  set -o errexit
}

get-chen-windows () { 
  cat ${chen_mchale_windows} | tail -n +2
}

get-clinvar-snvs () {
  cat ${clinvar_snvs} \
    | tail -n +2 \
    | cut -f1-3 \
    | sort --version-sort -k1,1 -k2,2n
}

add-clinvar-snv-counts () {
  bedtools intersect \
      -a <(get-chen-windows) \
      -b <(get-clinvar-snvs) \
      -c
}

(
  header-line 
  add-clinvar-snv-counts
) > ${chen_mchale_windows_with_clinvar_snv_counts}  

info "Wrote windows with clinvar snv counts to:" ${chen_mchale_windows_with_clinvar_snv_counts}  

################################ SPOT CHECKING #######################################

get-clinvar-snvs-full-records () {
  cat ${clinvar_snvs} \
    | tail -n +2 \
    | sort --version-sort -k1,1 -k2,2n
}

show-clinvar-variants-in-window () {
  window=$1
  info "clinvar SNVs in window:" "${window}"
  bedtools intersect \
      -a <(get-clinvar-snvs-full-records) \
      -b <(echo -e "${window}") 
}

show-clinvar-variants-in-window "chr11\t5227000\t5228000"


