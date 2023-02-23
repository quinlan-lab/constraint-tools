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
                                               
exon_adjacent_windows="${CONSTRAINT_TOOLS_DATA}/genome-wide-predictions/predict-germline-grch38.exon-adjacent-windows.bed.gz"
clinvar_snvs="${CONSTRAINT_TOOLS_DATA}/clinvar-noncoding-with-negative-controls/Clinvar_nc_snv_pathogenic.hg38.bed"
exon_adjacent_windows_with_clinvar_snv_counts="${CONSTRAINT_TOOLS_DATA}/genome-wide-predictions/predict-germline-grch38.exon-adjacent-windows-clinvar.bed"

header-line () {
  set +o errexit
  echo -e "$(zcat ${exon_adjacent_windows} | head -1)\tClinVar SNV count"
  set -o errexit
}

get-exon-adjacent-windows () { 
  zcat ${exon_adjacent_windows} | tail -n +2
}

get-clinvar-snvs () {
  cat ${clinvar_snvs} \
    | tail -n +2 \
    | cut -f1-3 \
    | sort --version-sort -k1,1 -k2,2n
}

add-clinvar-snv-counts () {
  bedtools intersect \
      -a <(get-exon-adjacent-windows) \
      -b <(get-clinvar-snvs) \
      -c
}

(
  header-line 
  add-clinvar-snv-counts
) > ${exon_adjacent_windows_with_clinvar_snv_counts}  

info "Wrote exon-adjacent windows with clinvar snv counts to:" ${exon_adjacent_windows_with_clinvar_snv_counts}  



