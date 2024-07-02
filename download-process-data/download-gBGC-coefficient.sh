set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset

source download-process-data/set-environment-variables.sh 

# no need to export PATH since it is already in the environment: 
# `printenv | grep -w PATH` returns non-zero output 
PATH="${CONSTRAINT_TOOLS}/bin:$PATH" 
PATH="${CONSTRAINT_TOOLS}/utilities:$PATH" 

# exporting these variables make them visible in python scripts below 
export GBGC_PATH="${CONSTRAINT_TOOLS_DATA}/GC-biased-gene-conversion"
export POP="EUR"

mkdir --parents ${GBGC_PATH}
info "download directory is:" ${GBGC_PATH}

download () { 
  # https://genome.cshlp.org/content/25/8/1215/suppl/DC1
  # https://genome.cshlp.org/content/suppl/2015/06/03/gr.185488.114.DC1/ReadMe.txt
  # Text S5 shows that coordinate system is hg18: 
  # "As a control, we also used the pedigree-based genetic map from deCODE (Kong et al. 2010). 
  #  As these maps are reported in the version hg18 of the human genome assembly, 
  #  we converted the location of SNPs from hg19 to hg18 coordinates using the liftOver tool 
  #  (https://genome.ucsc.edu/cgi-bin/hgLiftOver). 
  #  A small fraction of SNPs (N=31,056) could not be mapped onto hg18 and the SNPs were discarded."
  curl "https://genome.cshlp.org/content/suppl/2015/06/03/gr.185488.114.DC1/B_estimates_${POP}.txt" \
    | tr '\r' '\n' \
    | process_gBGC_coefficient_data \
    > ${GBGC_PATH}/gBGC-coefficient.hg18.${POP}.tsv
}

lift () {
  info "Lifting over..."
  bash \
    ${CONSTRAINT_TOOLS}/download-process-data/lift.sh \
    "${GBGC_PATH}/gBGC-coefficient.hg18.${POP}.bed" \
    hg18 \
    hg38 \
    "-bedPlus=3 -tab"
  mv \
    "${GBGC_PATH}/gBGC-coefficient.hg18.${POP}.bed.hg38" \
    "${GBGC_PATH}/gBGC-coefficient.hg38.${POP}.bed"
  mv \
    "${GBGC_PATH}/gBGC-coefficient.hg18.${POP}.bed.hg38.unmapped" \
    "${GBGC_PATH}/gBGC-coefficient.hg18.${POP}.umapped.bed"
}

add_header_line () {
  print_header_line \
    | cat - "${GBGC_PATH}/gBGC-coefficient.hg38.${POP}.bed" \
    > "${GBGC_PATH}/gBGC-coefficient.hg38.${POP}.header.bed"
}

# download
# prepare_for_liftover 
# lift 
add_header_line
