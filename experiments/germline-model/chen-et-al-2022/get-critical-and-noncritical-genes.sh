set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset 

source download-process-data/set-environment-variables.sh 

genes_path="${CONSTRAINT_TOOLS_DATA}/genes/grch38"

cat \
    ${genes_path}/autosomal-dominant.genes.sorted.bed \
    ${genes_path}/haploinsufficient.genes.sorted.bed \
  | sort -k1,1V -k2,2n \
  | uniq \
  > ${genes_path}/critical.genes.sorted.bed

get-noncritical-genes \
    ${genes_path}/critical.genes.sorted.bed \
    <(bgzip --decompress --stdout  ${genes_path}/genes.sorted.bed.gz) \
  > ${genes_path}/noncritical.genes.sorted.bed

info "Wrote noncritical genes to ${genes_path}/noncritical.genes.sorted.bed"

