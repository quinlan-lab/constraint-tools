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

chen_mchale_windows="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/mchale.kmerSizes.trainSets.noisy.enhancer-exon-cpgIsland.bed"

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5467550/
# data courtesy: Tom Nicholas 
GeneHancer_enhancers="/scratch/ucgd/lustre-work/quinlan/u0055382/genome_reference/genehancer/genehancer_GRCh38.bed.gz" 

get-chen-windows () { 
  cat ${chen_mchale_windows} | tail -n +2
}

get-GeneHancer-enhancers () {
  zcat ${GeneHancer_enhancers} \
    | tail -n +2 \
    | cut -f1-7 \
    | awk '{ print "chr"$0}' \
    | sort --version-sort -k1,1 -k2,2n \
    | uniq 
}

# get-chen-windows | head 
# get-GeneHancer-enhancers | head

header-line () {
  set +o errexit
  echo -e "$(head -1 ${chen_mchale_windows})\tenhancer_chromosome\tenhancer_start\tenhancer_end\tGHid\tenhancer_type\telite_enhancer\tgene_targeted_by_enhancer|enhancer_gene_association_score|elite_enhancer_gene_association\twindow_enhancer_overlap_bps"
  set -o errexit
}

add-enhancer-info () {
  bedtools intersect \
      -a <(get-chen-windows) \
      -b <(get-GeneHancer-enhancers) \
      -wao 
}

out_filename="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/mchale.kmerSizes.trainSets.noisy.overlapAmounts.cpg-islands.enhancer-info.bed"
(
  header-line
  add-enhancer-info
) > ${out_filename}
# | head -5 | column -t -s $'\t' 

info "Wrote:" ${out_filename}  

