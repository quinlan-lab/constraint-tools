
set -o errexit
set -o pipefail

set -o xtrace
# Must use single quote to prevent variable expansion.
# For example, if double quotes were used, ${LINENO} would take on the value of the current line,
# instead of its value when PS4 is used later in the script
# https://stackoverflow.com/a/6697845/6674256
# ${FOO:+val}    val if $FOO is set
# ${FOO[0]}   element #0 of the FOO array
# https://www.gnu.org/software/bash/manual/html_node/Bash-Variables.html
PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

export RED='\033[0;31m'
export CYAN='\033[0;36m'
export NO_COLOR='\033[0m' 

root="/scratch/ucgd/lustre-work/quinlan/u1240855/constraint-tools"
path="${root}/data/neutral_region"
gtf="${root}/data/genes/Homo_sapiens.GRCh37.74.sorted.gtf.gz"
gap="${root}/data/gap/gap.txt"

mkdir --parents ${path}

url="ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes"

download_file () {
    local url=$1
    wget ${url} --output-document=${path}/hg19.chrom.sizes
}

echo -e "${CYAN}Downloading UCSC hg19 chr size file...${NO_COLOR}\n"
download_file ${url}

echo -e "${CYAN}Selecting for chr1-22,X,Y and sorting...${NO_COLOR}\n"
for chr in {1..22} X Y; do grep -w "chr${chr}" $path/hg19.chrom.sizes ; done | sort -k1,1 > $path/hg19.chrom.sizes.chr.sorted

echo -e "${CYAN}Selecting for chr1-22,X,Y from the gtf file and sorting...${NO_COLOR}\n"
zcat $gtf | cut -f 1,4-5 | awk '{print "chr"$0""}' > $path/Homo_sapiens.GRCh37.74.sorted.bed
for chr in {1..22} X Y; do grep -w "chr${chr}" $path/Homo_sapiens.GRCh37.74.sorted.bed ; done | sort -k1,1 -k2,2n > $path/Homo_sapiens.GRCh37.74.chr.sorted.bed

echo -e "${CYAN}Sorting gap bed file...${NO_COLOR}\n"
cat $gap | cut -f 2-4,8 | sort -k1,1 -k2,2n > $path/gap_sorted.bed

echo -e "${CYAN}Generate bed file of regions putatively under neutral selection...${NO_COLOR}\n"
bedtools complement -i $path/Homo_sapiens.GRCh37.74.chr.sorted.bed -g $path/hg19.chrom.sizes.chr.sorted | bedtools intersect -a - -b $path/gap_sorted.bed -v > $path/putative_neutral_regions.bed

rm $$path/Homo_sapiens.GRCh37.74.sorted.bed
rm $path/hg19.chrom.sizes
