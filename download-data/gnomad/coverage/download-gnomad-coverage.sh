#!/bin/sh
#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw

## Define input variables
while [[ "$1" =~ ^- ]]; do
  case $1 in
    --coverage-threshold ) shift; [[ ! $1 =~ ^- ]] && coverage_threshold=$1;;
    --fraction-threshold ) shift; [[ ! $1 =~ ^- ]] && fraction_threshold=$1;;
    --version ) shift; [[ ! $1 =~ ^- ]] && version=$1;;
    --seq ) shift; [[ ! $1 =~ ^- ]] && seq=$1;;
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

info "Filterting on fraction-threshold value of: ${fraction_threshold}..."
info "Filtering for sites with 10x coverage..."
info "Filtering gnomad-${version} ${seq} coverage file to select for positions in which (${fraction_threshold} * 100)% of samples have a mean depth of value: 10x"

#######################################

if [ ${version} == "v3" ] && [ ${seq} == "wgs" ]; then
	url="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.0.1/coverage/genomes/gnomad.genomes.r3.0.1.coverage.summary.tsv.bgz"
	gnomad_coverage_file="gnomad_v3_coverage.summary.tsv.bgz"
	gnomad_coverage_filtered="gnomad_v3_coverage.filtered"

	## Define directory to download files into 
	coverage_path="${CONSTRAINT_TOOLS}/data/gnomad/v3/coverage"
	mkdir --parents ${coverage_path}	

	info "Downloading gnomad v3 wgs coverage file"
	#wget ${url} --output-document=${coverage_path}/${gnomad_coverage_file}
	#tail ${coverage_path}/${gnomad_coverage_file}	

	info "Applying above filter on gnomad v3 wgs coverage file..."
	## coverage column 7 for over 10x --> coverage_threshold=7
	zcat ${coverage_path}/${gnomad_coverage_file} | tail -n +2 | awk -v c=7 -v t=${fraction_threshold} '$c>t {print $1}' | sed 's/:/\t/g' | awk '{print $1"\t"($2-1)"\t"$2}' | bedtools merge > ${coverage_path}/${gnomad_coverage_filtered}.hg38.bed

## Download gnomad coverage data (v2)
elif [ ${seq} == "wgs" ] && [ ${version} == "v2" ]; then
	url="https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1/coverage/genomes/gnomad.genomes.coverage.summary.tsv.bgz"
	gnomad_coverage_file="gnomad_v2_wgs_coverage.summary.tsv.bgz"
	gnomad_coverage_filtered="gnomad_v2_wgs_coverage.filtered"
	
	## Define directory to download files into 
	coverage_path="${CONSTRAINT_TOOLS}/data/gnomad/v2/wgs/coverage"
	mkdir --parents ${coverage_path}

	info "Downloading gnomad v2 wgs coverage file"
	#wget ${url} --output-document=${coverage_path}/${gnomad_coverage_file}

	info "Performing the above filter..."
        ## coverage column 7 for over 10x --> coverage_threshold=7
	zcat ${coverage_path}/${gnomad_coverage_file} | tail -n+2 | awk -v c=7 -v t=${fraction_threshold} '$c>t {print "chr"$1"\t"($2-1)"\t"$2}' | bedtools merge > ${coverage_path}/${gnomad_coverage_filtered}.hg19.bed

	## https://gist.github.com/brentp/894555#file-lift-sh
	info "Perform liftover to convert hg19 coordinates to hg38"
	bash ${CONSTRAINT_TOOLS}/download-data/lift.sh ${coverage_path}/${gnomad_coverage_filtered}.hg19.bed hg19 hg38 

	## Remove intermediate files
	mv ${coverage_path}/${gnomad_coverage_filtered}.hg19.bed.hg38 ${coverage_path}/${gnomad_coverage_filtered}.hg38.bed
	rm ${coverage_path}/${gnomad_coverage_filtered}.hg19.bed

elif [ ${seq} == "wes" ] && [ ${version} == "v2" ]; then
        url="https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1/coverage/exomes/gnomad.exomes.coverage.summary.tsv.bgz"
        gnomad_coverage_file="gnomad_v2_wes_coverage.summary.tsv.bgz"
        gnomad_coverage_filtered="gnomad_v2_wes_coverage.filtered"

	## Define directory to download files into
	coverage_path="${CONSTRAINT_TOOLS}/data/gnomad/v2/wes/coverage"
	mkdir --parents ${coverage_path}
	
        info "Downloading gnomad v2 wes coverage file"
        #wget ${url} --output-document=${coverage_path}/${gnomad_coverage_file}
	
	info "Performing the above filter..."
        ## coverage column 7 for over 10x --> coverage_threshold=7
        zcat ${coverage_path}/${gnomad_coverage_file} | tail -n+2 | awk -v c=7 -v t=${fraction_threshold} '$c>t {print "chr"$1"\t"($2-1)"\t"$2}' | bedtools merge > ${coverage_path}/${gnomad_coverage_filtered}.hg19.bed	

	## https://gist.github.com/brentp/894555#file-lift-sh
	info "Perform liftover to convert hg19 coordinates to hg38"
        bash ${CONSTRAINT_TOOLS}/download-data/lift.sh ${coverage_path}/${gnomad_coverage_filtered}.hg19.bed hg19 hg38 

        ## Remove intermediate files
        mv ${coverage_path}/${gnomad_coverage_filtered}.hg19.bed.hg38 ${coverage_path}/${gnomad_coverage_filtered}.hg38.bed
        rm ${coverage_path}/${gnomad_coverage_filtered}.hg19.bed
else 
	echo "input valid gnomad version number and sequencing strategy employed..."
	exit 1
fi

info "Sorting and compressing coverage bed file..."
cat ${coverage_path}/${gnomad_coverage_filtered}.hg38.bed \
  | get-regular-chromosomes \
  | sort-compress-index-bed --name ${coverage_path}/${gnomad_coverage_filtered}.sorted


