set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset

source /scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools/set-environment-variables.sh 

lift () {
  local filename_root=$1
  bash \
    ${CONSTRAINT_TOOLS}/download-process-data/lift.sh \
    ${CONSTRAINT_TOOLS_DATA}/khurana/${filename_root}.hg19.sorted.coordinates-only.bed \
    hg19 \
    hg38
  mv \
    ${CONSTRAINT_TOOLS_DATA}/khurana/${filename_root}.hg19.sorted.coordinates-only.bed.hg38 \
    ${CONSTRAINT_TOOLS_DATA}/khurana/${filename_root}.hg38.hg19.bed

  mv \
    ${CONSTRAINT_TOOLS_DATA}/khurana/${filename_root}.hg19.sorted.coordinates-only.bed.hg38.unmapped \
    ${CONSTRAINT_TOOLS_DATA}/khurana/${filename_root}.hg19.unmapped.bed
}

lift "all-enhancers-with-network-features"
