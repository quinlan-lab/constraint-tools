set -o errexit
set -o pipefail
# set -o noclobber

set -o xtrace

source set-environment-variables.sh 

# p34 of: 
# https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-1969-6/MediaObjects/41586_2020_1969_MOESM1_ESM.pdf 
# says “used Ensembl release 74 as a base for the annotation”
genes_url="http://ftp.ensembl.org/pub/release-74/gtf/homo_sapiens"
genes="Homo_sapiens.GRCh37.74"

genes_path="${CONSTRAINT_TOOLS}/data/genes"

mkdir --parents ${genes_path}

info "Downloading GRCH37 annotations..."
wget ${genes_url}/${genes}.gtf.gz --output-document=${genes_path}/${genes}.gtf.gz

info "Sorting and block compressing GRCH37 annotations..."
set +o errexit
(
  zgrep "^#" ${genes_path}/${genes}.gtf.gz
  zgrep -v "^#" ${genes_path}/${genes}.gtf.gz | sort --version-sort -k1,1 -k4,4n
) | 
  bgzip > ${genes_path}/${genes}.sorted.gtf.gz
set -o errexit

info "Indexing GRCH37 annotations..." 
tabix ${genes_path}/${genes}.sorted.gtf.gz

