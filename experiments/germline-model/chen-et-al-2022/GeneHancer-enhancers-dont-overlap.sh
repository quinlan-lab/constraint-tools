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

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5467550/
# data courtesy: Tom Nicholas 
GeneHancer_enhancers="/scratch/ucgd/lustre-work/quinlan/u0055382/genome_reference/genehancer/genehancer_GRCh38.bed.gz" 

get-GeneHancer-enhancers () {
  zcat ${GeneHancer_enhancers} \
    | tail -n +2 \
    | cut -f1-3 \
    | awk -v OFS="\t" '{ print "chr"$1, $2, $3 }' \
    | sort --version-sort -k1,1 -k2,2n \
    | uniq 
}

number_enhancers=$(get-GeneHancer-enhancers | wc -l)

# get-GeneHancer-enhancers | head 

find-overlaps () {
  bedtools intersect \
      -a <(get-GeneHancer-enhancers) \
      -b <(get-GeneHancer-enhancers) \
      -wa -wb
}

number_intersections=$(find-overlaps | wc -l)

if [[ ${number_intersections} > ${number_enhancers} ]]; then 
  info "There are more intersections than enhancers!"
else
  info "Enhancers are non-overlapping" 
fi 


