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

WINDOW_SIZE="$1" 
FEATURES_FILENAME="$2"
COVERAGE_FILENAME="$3"
SLIDING_WINDOW_SIZE="$4"

step_size=$((SLIDING_WINDOW_SIZE/2))

make_sliding_windows () {
  bedtools makewindows \
    -g <(echo -e "synthetic_chrom\t${WINDOW_SIZE}") \
    -w ${SLIDING_WINDOW_SIZE} \
    -s ${step_size} 
}

bedtools coverage \
  -a <(make_sliding_windows) \
  -b $FEATURES_FILENAME \
  -counts \
  > ${COVERAGE_FILENAME}

info "Wrote coverage:" "${COVERAGE_FILENAME}"
