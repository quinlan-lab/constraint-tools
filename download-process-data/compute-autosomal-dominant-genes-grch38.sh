set -o errexit
set -o pipefail
# set -o noclobber

# set -o xtrace
set -o nounset 

source download-process-data/set-environment-variables.sh 

genes_path="${CONSTRAINT_TOOLS_DATA}/genes/grch38"

# https://github.com/macarthur-lab/gene_lists
autosomal_dominant_gene_symbols="https://raw.githubusercontent.com/macarthur-lab/gene_lists/master/lists/berg_ad.tsv" 

info "Merging the autosomal dominant genes with all genes to extract the coordinates of the automosomal dominant genes..." 
get-gene-coordinates-from-symbols \
    <(curl ${autosomal_dominant_gene_symbols}) \
    <(zcat ${genes_path}/genes.sorted.bed.gz) \
  | get-regular-chromosomes \
  | sort -k1,1V -k2,2n \
> ${genes_path}/autosomal-dominant.genes.sorted.bed
info "Wrote" "${genes_path}/autosomal-dominant.genes.sorted.bed"


