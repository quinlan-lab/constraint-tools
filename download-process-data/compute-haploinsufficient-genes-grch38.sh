set -o errexit
set -o pipefail
# set -o noclobber

# set -o xtrace
set -o nounset 

source download-process-data/set-environment-variables.sh 

genes_path="${CONSTRAINT_TOOLS_DATA}/genes/grch38"

# https://github.com/macarthur-lab/gene_lists
haploinsufficient_gene_symbols="https://raw.githubusercontent.com/macarthur-lab/gene_lists/master/lists/clingen_level3_genes_2018_09_13.tsv" 

info "Merging the haploinsufficient genes with all genes to extract the coordinates of the haploinsufficient genes..." 
get-gene-coordinates-from-symbols \
    <(curl ${haploinsufficient_gene_symbols}) \
    <(zcat ${genes_path}/genes.sorted.bed.gz) \
  | get-regular-chromosomes \
  | sort -k1,1V -k2,2n \
> ${genes_path}/haploinsufficient.genes.sorted.bed
info "Wrote" "${genes_path}/haploinsufficient.genes.sorted.bed"


