set -o errexit
set -o pipefail
# set -o xtrace
set -o nounset

source download-process-data/set-environment-variables.sh 

mendelian_variants="${CONSTRAINT_TOOLS_DATA}/noncoding-variants-associated-with-Mendelian-traits/noncoding-mendelian-variants.hg38.txt"
canonical_exons="${CONSTRAINT_TOOLS_DATA}/genes/grch38/canonical-exons.sorted.bed.gz"

get-mendelian-variants-full-record () {
  cat ${mendelian_variants} \
    | tail -n +2 \
    | python ${CONSTRAINT_TOOLS}/experiments/germline-model/chen-et-al-2022/convert_to_bed.py \
    | sort --version-sort -k1,1 -k2,2n
}

show-mendelian-variants-in-canonical-exons () {
  bedtools intersect \
      -a <(get-mendelian-variants-full-record) \
      -b ${canonical_exons} \
      -wa -wb
}

info "Mendelian variants in canonical exons:" 
show-mendelian-variants-in-canonical-exons

info "number of variants:" $(get-mendelian-variants-full-record | wc -l)
info "number of variants that overlap canonical exons:" $(show-mendelian-variants-in-canonical-exons | tail -n +2 | wc -l)
