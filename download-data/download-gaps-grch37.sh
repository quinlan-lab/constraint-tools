set -o errexit
set -o pipefail
# set -o noclobber
set -o xtrace
set -o nounset

source download-data/set-environment-variables.sh 

# description & schema: 
# https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1139554535_L3a5g9m3mEN7R4zFxGxarD9lWCcB&jsh_pageVertPos=0&clade=mammal&org=Human&db=hg19&hgta_group=map&hgta_track=gap&hgta_table=gap&hgta_doSchema=describe+table+schema&hgta_regionType=genome&position=chrX%3A15%2C578%2C261-15%2C621%2C068&hgta_outputType=primaryTable&boolshad.sendToGalaxy=0&boolshad.sendToGreat=0&hgta_outFileName=&hgta_compressType=none
gaps_url="ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/"
gaps_data="gap.txt.gz"
gaps_schema="gap.sql"

# data now stored at /scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/gaps/grch37/
gaps_path="${CONSTRAINT_TOOLS}/data/gaps"

mkdir --parents ${gaps_path}

info "Downloading gaps..."
wget ${gaps_url}/${gaps_schema} --output-document=${gaps_path}/${gaps_schema}
curl --location ${gaps_url}/${gaps_data} \
  | cut -f 2-4 \
  | get-regular-chromosomes \
  | sort-compress-index-bed --name ${gaps_path}/gaps.sorted



