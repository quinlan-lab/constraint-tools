#!/bin/sh
#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --nodes=3
#SBATCH --ntasks=16
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH -o logs/download-gnomad-coverage.out

## Define input variables
while [[ "$1" =~ ^- ]]; do
  case $1 in
    --coverage ) shift; [[ ! $1 =~ ^- ]] && coverage=$1;;
    --threshold ) shift;  [[ ! $1 =~ ^- ]] && threshold=$1;;
    *) error "$0: $1 is an invalid flag"; exit 1;;
  esac
  shift
done

set -o errexit
set -o pipefail
# set -o noclobber
set -o xtrace
set -o nounset

source download-data/set-environment-variables.sh 

## Define directory to download files into 
gnomad_path="${CONSTRAINT_TOOLS}/data/gnomad/v3/coverage"

mkdir --parents ${gnomad_path}
 
echo "Threshold: ${threshold}"
 
## Specify threshold
if [ ${threshold} -gt 100 ]; then	
	echo "Detected threshold value greater than 100... Please supply value that is a fraction..."
	exit 1

elif [ ${threshold} -lt 1 ]; then
	echo "Detected threshold value less than 1... Using this value as a fraction (e.g. 0.8 --> 80%)..."

elif [ ${threshold} -gt 1 ]; then
	echo "Detected threshold value greater than 1 --> attempting to convert to a fraction (e.g. 80% --> 0.8)..."
	threshold=$(echo "${threshold} / 100" | bc -l)

else 
	echo "Invalid threshold value provided..." 
	exit 1
fi
 
## Specicy column to filter on based on desired coverage parameter
if [ ${coverage} == 1 ]; then
	column=5

elif [ ${coverage} == 5 ]; then
        column=6

elif [ ${coverage} == 10 ]; then
        column=7

elif [ ${coverage} == 15 ]; then
        column=8

elif [ ${coverage} == 20 ]; then
	column=9

elif [ ${coverage} == 25 ]; then
	column=10

elif [ ${coverage} == 30 ]; then
	column=11

elif [ ${coverage} == 50 ]; then
	column=12

elif [ ${coverage} == 100 ]; then
	column=13

else
	echo "input valid coverage value (1, 5, 10, 15, 20, 25, 50, or 100)..."
	exit 1
fi

info "Filtering gnomad v3 coverage file to select for positions in which (${threshold} * 100)% of samples have a mean depth of value: ${coverage}X"

## Download gnomad coverage data (v3)
url="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.0.1/coverage/genomes/gnomad.genomes.r3.0.1.coverage.summary.tsv.bgz"
gnomad_coverage="gnomad_v3_coverage.summary.tsv.bgz"
gnomad_coverage_filtered="gnomad_v3_coverage_${coverage}X_${threshold}"

info "Downloading gnomad v3 coverage file"
wget ${url} --output-document=${gnomad_path}/${gnomad_coverage}

info "Performing the above filter..."
zcat ${gnomad_path}/${gnomad_coverage} --force | tail -n+2 | awk -v c=${coverage} -v t=${threshold} '$c>t {print $1}' | sed 's/:/\t/g' | awk '{print $1"\t"($2-1)"\t"$2}' | bedtools merge > ${gnomad_path}/${gnomad_coverage_filtered}.bed

info "Sorting and compressing coverage bed file..."
cat ${gnomad_coverage_filtered} | sort-compress-index-bed --name ${gnomad_coverage_filtered}.sorted.bed.gz
