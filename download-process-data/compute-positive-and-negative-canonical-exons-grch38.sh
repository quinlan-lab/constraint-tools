set -o errexit
set -o pipefail
# set -o noclobber

# set -o xtrace
set -o nounset 

source download-process-data/set-environment-variables.sh 

genes_path="${CONSTRAINT_TOOLS_DATA}/genes/grch38"

# https://github.com/macarthur-lab/gene_lists
autosomal_dominant_gene_symbols="https://raw.githubusercontent.com/macarthur-lab/gene_lists/master/lists/berg_ad.tsv" 
haploinsufficient_gene_symbols="https://raw.githubusercontent.com/macarthur-lab/gene_lists/master/lists/clingen_level3_genes_2018_09_13.tsv"
olfactor_receptor_gene_symbols="https://raw.githubusercontent.com/macarthur-lab/gene_lists/master/lists/olfactory_receptors.tsv"

get-header () {
  local header="chromosome\texon_start\texon_end\tgene_symbol\texon_rank\tgene_biotype"
  echo -e "${header}" 
}

get-positive-gene-symbols () { 
  (
    curl ${autosomal_dominant_gene_symbols}
    curl ${haploinsufficient_gene_symbols}
  ) | sort | uniq
}

get-negative-gene-symbols () { 
  curl ${olfactor_receptor_gene_symbols} | sort | uniq
}

write-canonical-exons-in-class () { 
  local class=$1

  local output="${genes_path}/canonical-exons.${class}.sorted.bed"
  info "Merging the gene symbols of the ${class} genes with all canonical exons to extract the coordinates of the ${class} canonical exons..." 
  (
    get-header 
    get-exon-coordinates-from-symbols \
        <(get-${class}-gene-symbols) \
        <(zcat ${genes_path}/canonical-exons.sorted.bed.gz) \
        <(get-header) \
      | get-regular-chromosomes \
      | sort -k1,1V -k2,2n 
  ) > ${output}
  info "Wrote" ${output}
}

write-canonical-exons-in-class "positive"
write-canonical-exons-in-class "negative"
