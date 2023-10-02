set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset

source download-process-data/set-environment-variables.sh 

# schema:
# https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1720914064_N7rJbaDuFwHPsvFPPk141LzwBzeA&clade=mammal&org=Human&db=hg38&hgta_group=allTracks&hgta_track=cytoBand&hgta_table=0&hgta_regionType=genome&position=chr11%3A134%2C074%2C275-134%2C080%2C274&hgta_outputType=primaryTable&hgta_outFileName=

url="ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/database"        
data="cytoBand.txt.gz"

dir="${CONSTRAINT_TOOLS_DATA}/chromosome-bands/grch38"

mkdir --parents ${dir}

info "Downloading chromosome bands..."
curl --location ${url}/${data} \
  | gunzip -c \
  | get-regular-chromosomes \
  | get-nonXY-chromosomes \
  | sort -k1,1 -k2,2n --version-sort \
  > ${dir}/chromosome-bands.sorted.bed

info "Find and store centromere locations..." 
grep --ignore-case "cen" ${dir}/chromosome-bands.sorted.bed \
  | awk '{print $1,$2,$3}' \
  | sort -k1,1 -k2,2n --version-sort \
  > ${dir}/centromeres.bed  

info "Find and store telomere locations (p-arm and q-arm)..."
awk '{print $1,$2,$3}' ${dir}/chromosome-bands.sorted.bed \
  | extract-telomeres \
  | sort -k1,1 -k2,2n --version-sort \
  > ${dir}/telomeres.bed

