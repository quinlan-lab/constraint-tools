set -o errexit
set -o pipefail
# set -o noclobber

set -o xtrace
set -o nounset 

source download-data/set-environment-variables.sh 

# gff3 format spec: 
# https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
genes_url="http://ftp.ensembl.org/pub/current_gff3/homo_sapiens/"           
genes="Homo_sapiens.GRCh38.105"
README="README"

genes_path="${CONSTRAINT_TOOLS_DATA}/genes/grch38"
mkdir --parents ${genes_path}

#######################

# info "Downloading GRCH38 annotations..."
# wget ${genes_url}/${genes}.gff3.gz --output-document=${genes_path}/${genes}.gff3.gz
# wget ${genes_url}/${README} --output-document=${genes_path}/${README}

info "Extracting and processing exons..."
zcat ${genes_path}/${genes}.gff3.gz \
  | grep -v "^#" \
  | awk --assign OFS='\t' '$3 == "exon" { print "chr"$1, $4, $5 }' \
  | get-regular-chromosomes \
  | sort-compress-index-bed --name ${genes_path}/exons.sorted


