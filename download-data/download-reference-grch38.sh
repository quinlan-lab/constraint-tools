#!/usr/bin/env/bash
set -o errexit
set -o pipefail
# set -o noclobber
set -o xtrace
set -o nounset 

source set-environment-variables.sh 

# https://www.gungorbudak.com/blog/2018/05/16/how-to-download-hg38-grch38-fasta-human-reference-genome/
# http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use

reference_url="http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet"
reference_genome="hg38.analysisSet"

# data now stored at: /scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/reference/grch38
reference_path=${CONSTRAINT_TOOLS}/data/reference/grch38

mkdir --parents ${reference_path}

download_file_with_suffix () {
	local suffix_=$1
	wget ${reference_url}/${reference_genome}.${suffix_} --output-document=${reference_path}/${reference_genome}.${suffix_}
}

info "Downloading GRCh38 reference that excludes ALT contigs..."
wget ${reference_url}/README.txt --output-document=${reference_path}/README.txt
#download_file_with_suffix "fa.gz"

check_digest () {
	local suffix_=$1
	local expected_digest_=$2
	observed_digest_=$(md5sum ${reference_path}/${reference_genome}.${suffix_} | awk '{ print $1 }')

	if [[ ${observed_digest_} == ${expected_digest_} ]]; then
		info "${reference_path}/${reference_genome}.${suffix_}: digest check passed"
	else 
		info "${reference_path}/${reference_genome}.${suffix_}: digest check failed"
	fi
} 

check_digest "fa.gz" "813545f7c55ce7a9bb6c2de81321f242" ## check that this is the correct digest

info "Decompressing reference file..."
# https://unix.stackexchange.com/a/158538/406037
set +o errexit
gzip --decompress ${reference_path}/${reference_genome}.fa.gz 
set -o errexit

info "Block compressing reference file..."
bgzip --stdout ${reference_path}/${reference_genome}.fa > ${reference_path}/${reference_genome}.fa.gz

info "Indexing block-compressed fasta file...."
samtools faidx ${reference_path}/${reference_genome}.fa.gz

info "Removing uncompressed fasta file ... "
rm ${reference_path}/${reference_genome}.fa

