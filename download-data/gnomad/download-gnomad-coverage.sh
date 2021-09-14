#!/bin/sh
#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH -o logs/download-gnomad-coverage.out

## Define input variables
while [[ "$1" =~ ^- ]]; do
  case $1 in
    --coverage ) shift; [[ ! $1 =~ ^- ]] && coverage=$1;;
    --threshold ) shift; [[ ! $1 =~ ^- ]] && threshold=$1;;
    --version ) shift; [[ ! $1 =~ ^- ]] && version=$1;;
    --sequencing ) shift; [[ ! $1 =~ ^- ]] && sequencing=$1;;
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

info "Filtering gnomad-${version} ${sequencing} coverage file to select for positions in which (${threshold} * 100)% of samples have a mean depth of value: ${coverage}X"

#######################################

if [ ${version} == "v3" ] && [ ${sequencing} == "wgs" ]; then
	url="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.0.1/coverage/genomes/gnomad.genomes.r3.0.1.coverage.summary.tsv.bgz"
	gnomad_coverage_file="gnomad_v3_coverage.summary.tsv.bgz"
	gnomad_coverage_filtered="gnomad_v3_coverage.filtered"

	## Define directory to download files into 
	coverage_path="${CONSTRAINT_TOOLS}/data/gnomad/v3/coverage"
	mkdir --parents ${coverage_path}	

	info "Downloading gnomad v3 wgs coverage file"
	#wget ${url} --output-document=${coverage_path}/${gnomad_coverage_file}
	#tail ${coverage_path}/${gnomad_coverage_file}	

	info "Performing the above filter..."
	zless ${coverage_path}/${gnomad_coverage_file} | tail -n +2 | awk -v c=${column} -v t=${threshold} '$c>t {print $1}' | sed 's/:/\t/g' | awk '{print $1"\t"($2-1)"\t"$2}' | bedtools merge > ${coverage_path}/${gnomad_coverage_filtered}.bed.hg38

## Download gnomad coverage data (v2)
elif [ ${sequencing} == "wgs" ] && [ ${version} == "v2" ]; then
	url="https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1/coverage/genomes/gnomad.genomes.coverage.summary.tsv.bgz"
	gnomad_coverage_file="gnomad_v2_wgs_coverage.summary.tsv.bgz"
	gnomad_coverage_filtered="gnomad_v2_wgs_coverage.filtered"
	
	## Define directory to download files into 
	coverage_path="${CONSTRAINT_TOOLS}/data/gnomad/v2/wgs/coverage"
	mkdir --parents ${coverage_path}

	info "Downloading gnomad v2 wgs coverage file"
	wget ${url} --output-document=${coverage_path}/${gnomad_coverage_file}

	info "Performing the above filter..."
	zless ${coverage_path}/${gnomad_coverage_file} | tail -n+2 | awk -v c=${column} -v t=${threshold} '$c>t {print "chr"$1"\t"($2-1)"\t"$2}' | bedtools merge > ${coverage_path}/${gnomad_coverage_filtered}.bed

	## https://gist.github.com/brentp/894555#file-lift-sh
	info "Perform liftover to convert hg19 coordinates to hg38"
	bash ${CONSTRAINT_TOOLS}/download-data/lift.sh ${coverage_path}/${gnomad_coverage_filtered}.bed hg19 hg38

elif [ ${sequencing} == "wes" ] && [ ${version} == "v2" ]; then
        url="https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1/coverage/exomes/gnomad.exomes.coverage.summary.tsv.bgz"
        gnomad_coverage_file="gnomad_v2_wes_coverage.summary.tsv.bgz"
        gnomad_coverage_filtered="gnomad_v2_wes_coverage.filtered"

	## Define directory to download files into
	coverage_path="${CONSTRAINT_TOOLS}/data/gnomad/v2/wes/coverage"
	mkdir --parents ${coverage_path}
	
        info "Downloading gnomad v2 wes coverage file"
        wget ${url} --output-document=${coverage_path}/${gnomad_coverage_file}
	
	info "Performing the above filter..."
        zless ${coverage_path}/${gnomad_coverage_file} | tail -n+2 | awk -v c=${column} -v t=${threshold} '$c>t {print "chr"$1"\t"($2-1)"\t"$2}' | bedtools merge > ${coverage_path}/${gnomad_coverage_filtered}.bed
	
	## https://gist.github.com/brentp/894555#file-lift-sh
	info "Perform liftover to convert hg19 coordinates to hg38"
        bash ${CONSTRAINT_TOOLS}/download-data/lift.sh ${coverage_path}/${gnomad_coverage_filtered}.bed hg19 hg38
else 
	echo "input valid gnomad version number and sequencing strategy employed..."
	exit 1
fi

info "Sorting and compressing coverage bed file..."
cat ${coverage_path}/${gnomad_coverage_filtered}.bed.hg38 \
  | get-regular-chromosomes \
  | sort-compress-index-bed --name ${coverage_path}/${gnomad_coverage_filtered}.sorted







