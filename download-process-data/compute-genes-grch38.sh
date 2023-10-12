set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset 

source download-process-data/set-environment-variables.sh 

genes_path="${CONSTRAINT_TOOLS_DATA}/genes/grch38"

# c.f. download-process-data/download-genes-grch38.sh
annotations="${genes_path}/Homo_sapiens.GRCh38.105.gff3.gz"

info "Extracting the coordinates and symbols of genes..."
zcat ${annotations} \
  | get-gene-coordinates-and-symbols \
  | get-regular-chromosomes \
  | sort-compress-index-bed --name ${genes_path}/genes.sorted

