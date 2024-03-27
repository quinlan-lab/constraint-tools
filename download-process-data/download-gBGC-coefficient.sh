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
export gBGC_path="${CONSTRAINT_TOOLS_DATA}/GC-biased-gene-conversion"
export pop="EUR"

mkdir --parents ${gBGC_path}
info "download directory is:" ${gBGC_path}

# https://genome.cshlp.org/content/25/8/1215/suppl/DC1
# https://genome.cshlp.org/content/suppl/2015/06/03/gr.185488.114.DC1/ReadMe.txt
curl "https://genome.cshlp.org/content/suppl/2015/06/03/gr.185488.114.DC1/B_estimates_${pop}.txt" \
  | tr '\r' '\n' \
  | process_gBGC_coefficient_data \
  > ${gBGC_path}/gBGC-coefficient.hg18.${pop}.tsv

prepare_for_liftover 

lift () {
  info "Lifting over..."
  bash \
    ${CONSTRAINT_TOOLS}/download-process-data/lift.sh \
    "${gBGC_path}/gBGC-coefficient.hg18.EUR.bed" \
    hg18 \
    hg38 \
    "-bedPlus=3 -tab"
  mv \
    "${gBGC_path}/gBGC-coefficient.hg18.EUR.bed.hg38" \
    "${gBGC_path}/gBGC-coefficient.hg38.EUR.bed"
  mv \
    "${gBGC_path}/gBGC-coefficient.hg18.EUR.bed.hg38.unmapped" \
    "${gBGC_path}/gBGC-coefficient.hg18.EUR.umapped.bed"
}

lift 

