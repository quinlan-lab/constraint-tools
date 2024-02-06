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

DATA_DIRECTORY="${1}"
DATA_STEM="${2}"
PUBLIC_REPO_DIR="${3}"
TRACK_NAME="${4}"
TRACK_DESCRIPTION="${5}"

create-bedgraph () {
  # https://genome.ucsc.edu/goldenPath/help/bedgraph.html
  (
    echo -e \
"track \
type=bedGraph \
name=${TRACK_NAME} \
description=${TRACK_DESCRIPTION} \
visibility=full \
color=0,0,0" 
    set +o errexit
    sort -k1,1V -k2,2n  "${DATA_DIRECTORY}/${DATA_STEM}.bed"
    set -o errexit
  ) > "${PUBLIC_REPO_DIR}/${DATA_STEM}.bedGraph" 
  info \
    "Wrote:" \
    "${PUBLIC_REPO_DIR}/${DATA_STEM}.bedGraph" 
}

push-to-public-repo () { 
  cd ${PUBLIC_REPO_DIR}
  git lfs track "${DATA_STEM}.bedGraph"
  git add "${DATA_STEM}.bedGraph" ".gitattributes"
  ( git commit -m "Add ${DATA_STEM} in UCSC-genome-browser format (bedgraph)" ) || true
  git push
  git lfs prune --recent
  info "Pushed ${DATA_STEM}.bedGraph to public repo"
}

create-bedgraph
push-to-public-repo
