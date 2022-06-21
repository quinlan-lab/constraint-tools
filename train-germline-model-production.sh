set -o errexit
set -o pipefail
set -o nounset

# sbatch train-germline-model.sh --constraint-tools-directory $PWD  

sbatch train-germline-model-production-window-sizes.sh 
