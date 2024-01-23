set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset

source /scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools/set-environment-variables.sh 

PATH="${CONSTRAINT_TOOLS}/experiments/germline-model/chen-et-al-2022:$PATH"

# https://stackoverflow.com/a/43476575/6674256
# need to export PYTHONPATH since it is not already in the environment: 
# `printenv | grep -w PYTHONPATH` returns zero output
export PYTHONPATH="${CONSTRAINT_TOOLS}/utilities"

CHEN_DATA_DIRECTORY="${1}"
CHEN_FILE_STEM="${2}"
PUBLIC_REPO_DIR="${3}"

create-bedgraph () {
  # https://genome.ucsc.edu/goldenPath/help/bedgraph.html
  (
    echo -e \
"track \
type=bedGraph \
name=chen-zscores \
description=published-chen-zscores-per-kb \
visibility=full \
color=0,0,255" 
    set +o errexit
    sort -k1,1V -k2,2n  "${CHEN_DATA_DIRECTORY}/${CHEN_FILE_STEM}.bed"
    set -o errexit
  ) > "${CHEN_DATA_DIRECTORY}/${CHEN_FILE_STEM}.bedGraph"
  info \
    "Wrote Chen z-scores in UCSC-genome-browser format to:" \
    "${CHEN_DATA_DIRECTORY}/${CHEN_FILE_STEM}.bedGraph"
}

push-to-public-repo () { 
  mv "${CHEN_DATA_DIRECTORY}/${CHEN_FILE_STEM}.bedGraph" "${PUBLIC_REPO_DIR}/${CHEN_FILE_STEM}.bedGraph" 
  cd ${PUBLIC_REPO_DIR}
  git add ${CHEN_FILE_STEM}.bedGraph
  ( git commit -m "Add Chen zscores in UCSC-genome-browser format (bedgraph)" ) || true
  git push
  info "Pushed ${CHEN_FILE_STEM}.bedGraph to public repo"
}

create-bedgraph
push-to-public-repo
