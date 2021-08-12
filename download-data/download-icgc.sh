set -o errexit
set -o pipefail
# set -o noclobber

set -o xtrace

source set-environment-variables.sh 

# https://dcc.icgc.org/releases/PCAWG/consensus_snv_indel/
# https://www.nature.com/articles/s41586-020-1969-6#MOESM3 : 
# 1. Supplementary Table 4
# 2. "Matching normal samples were obtained from blood (2,064 donors)" 
icgc_url="https://dcc.icgc.org/api/v1/download?fn=/PCAWG/consensus_snv_indel/final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz"
icgc_path="${CONSTRAINT_TOOLS}/data/icgc"

mkdir --parents ${icgc_path}

info "Downloading ICGC mutations..."
wget ${icgc_url} --output-document=${icgc_path}/mutations.maf.gz 

stream_file () {
  zgrep -v "^#" ${icgc_path}/mutations.maf.gz | 
    # remove carriage returns: 
    # https://stackoverflow.com/questions/800030/remove-carriage-return-in-unix
    dos2unix 
}

info "Sorting and block compressing ICGC mutations..."
set +o errexit
(
  head -1 <(stream_file)
  tail -n +2 <(stream_file) |
    sort --version-sort -k2,2 -k3,3n -k4,4n
) |
  bgzip > ${icgc_path}/mutations.sorted.maf.gz
set -o errexit

info "Indexing ICGC mutations..."
# http://www.htslib.org/doc/tabix.html
tabix \
    --skip-lines 1 \
    --sequence 2 \
    --begin 3 \
    --end 4 \
    --force \
  ${icgc_path}/mutations.sorted.maf.gz