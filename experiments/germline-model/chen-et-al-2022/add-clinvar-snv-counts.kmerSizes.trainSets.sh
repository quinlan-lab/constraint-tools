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
                                               
windows="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/chen-mchale.kmerSizes.trainSets.enhancer-exon.bed"
clinvar_snvs="${CONSTRAINT_TOOLS_DATA}/clinvar-noncoding-with-negative-controls/Clinvar_nc_snv_pathogenic.hg38.bed"
windows_with_clinvar_snv_counts="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/chen-mchale.kmerSizes.trainSets.enhancer-exon.clinvar.bed"

header-line () {
  set +o errexit
  echo -e "$(cat ${windows} | head -1)\tClinVar SNV count"
  set -o errexit
}

get-windows () { 
  cat ${windows} | tail -n +2
}

get-clinvar-snvs () {
  cat ${clinvar_snvs} \
    | tail -n +2 \
    | cut -f1-3 \
    | sort --version-sort -k1,1 -k2,2n
}

add-clinvar-snv-counts () {
  bedtools intersect \
      -a <(get-windows) \
      -b <(get-clinvar-snvs) \
      -c
}

(
  header-line 
  add-clinvar-snv-counts
) > ${windows_with_clinvar_snv_counts}  

info "Wrote windows with clinvar snv counts to:" ${windows_with_clinvar_snv_counts}  



