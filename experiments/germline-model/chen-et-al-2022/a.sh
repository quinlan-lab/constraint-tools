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
# https://docs.google.com/presentation/d/1o1rOcvnj69F9-qfXK0r19UKLupoivBvBlzBboM9u2dQ/edit#slide=id.g22aad12e963_0_24
GeneHancer_enhancers="/scratch/ucgd/lustre-work/quinlan/u0055382/genome_reference/genehancer/genehancer_GRCh38.bed.gz" 

# TODO: 
# create windows by subtracting 200bp and adding 400bp to enhancer_end (200bp < mean enhancer length (700bp))
get-windows-tail () {
  zcat ${GeneHancer_enhancers} \
    | tail -n +2 \
    | cut -f1-7 \
    | awk '{ print "chr"$0}' \
    | sort --version-sort -k1,1 -k2,2n \
    | uniq 
}

# TODO: 
# copy experiments/germline-model/chen-et-al-2022/intersect-windows-with-deletions.sh here 
# edit the new code to do: 
# 1. Bedtools-intersect 600bp windows (keep original enhancer coordinates for later use) and deletions (keep singleton-status of each deletion)
# 2. For each record, subtract center of window from deletion coordinates, and truncate the coordinates so that they lie in (-200, 400)
# 3. write altered deletion coordinates to disk 

