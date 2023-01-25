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
mendelian_variants="${CONSTRAINT_TOOLS_DATA}/noncoding-variants-associated-with-Mendelian-traits/noncoding-mendelian-variants.hg38.txt"
chen_mchale_windows_with_mendelian_variant_counts="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/chen-mchale-mendelian-variants.bed"

header-line () {
  set +o errexit
  echo -e "$(head -1 ${chen_mchale_windows})\tMendelian variant count"
  set -o errexit
}

get-chen-windows () { 
  cat ${chen_mchale_windows} | tail -n +2
}

get-mendelian-variants () {
  cat ${mendelian_variants} \
    | tail -n +2 \
    | cut -f1,2 \
    | awk -v OFS=$'\t' '{ print $1, $2-1, $2 }' \
    | sort --version-sort -k1,1 -k2,2n
}

add-mendelian-variant-counts () {
  bedtools intersect \
      -a <(get-chen-windows) \
      -b <(get-mendelian-variants) \
      -c
}

(
  header-line 
  add-mendelian-variant-counts
) > ${chen_mchale_windows_with_mendelian_variant_counts}  

info "Wrote windows with mendelian variant counts to:" ${chen_mchale_windows_with_mendelian_variant_counts}  


