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

windows="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/mchale.kmerSizes.trainSets.noisy.enhancer-exon.bed"
positive_vista_enhancers="${CONSTRAINT_TOOLS_DATA}/vista-enhancers/vista-enhancers.positive.hg38.hg19.tsv"
windows_with_positive_vista_enhancer_counts="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/mchale.kmerSizes.trainSets.noisy.enhancer-exon.positive-vista-enhancers.bed"

header-line () {
  set +o errexit
  echo -e "$(head -1 ${windows})\tpositive-vista-enhancer count"
  set -o errexit
}

get-windows-tail () { 
  cat ${windows} | tail -n +2
}

get-positive-vista-enhancers () {
  cut -f1-3 ${positive_vista_enhancers} \
    | sort --version-sort -k1,1 -k2,2n \
    | uniq 
}

add-positive-vista-enhancer-counts () {
  bedtools intersect \
      -a <(get-windows-tail) \
      -b <(get-positive-vista-enhancers) \
      -c
      # -wo 
}

(
  header-line 
  add-positive-vista-enhancer-counts
) > ${windows_with_positive_vista_enhancer_counts}  

info "Wrote positive vista enhancer counts to" ${windows_with_positive_vista_enhancer_counts}  


