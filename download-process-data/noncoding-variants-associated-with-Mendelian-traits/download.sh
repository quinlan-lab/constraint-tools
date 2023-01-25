set -o errexit
set -o pipefail
set -o xtrace
set -o nounset

source download-process-data/set-environment-variables.sh 

# https://www.nature.com/articles/s41588-018-0062-7
# "The noncoding variants associated with Mendelian traits (n = 427; used in Fig. 3b) are provided in Supplementary Table 3."
url='https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0062-7/MediaObjects/41588_2018_62_MOESM5_ESM.txt'
data_directory="/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/noncoding-variants-associated-with-Mendelian-traits"

curl ${url} > $data_directory/noncoding-mendelian-variants.hg38.txt

