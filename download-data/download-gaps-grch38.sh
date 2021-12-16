set -o errexit
set -o pipefail
# set -o noclobber
set -o xtrace
set -o nounset

source download-data/set-environment-variables.sh 

# https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1139554535_L3a5g9m3mEN7R4zFxGxarD9lWCcB&jsh_pageVertPos=0&clade=mammal&org=Human&db=hg38&hgta_group=map&hgta_track=gap&hgta_table=gap&hgta_doSchema=describe+table+schema
gaps_url="ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/database"        
gaps_data="gap.txt.gz"
gaps_schema="gap.sql"

gaps_path="${CONSTRAINT_TOOLS_DATA}/gaps/grch38"

mkdir --parents ${gaps_path}

info "Downloading gaps..."
wget ${gaps_url}/${gaps_schema} --output-document=${gaps_path}/${gaps_schema}
curl --location ${gaps_url}/${gaps_data} \
  | gunzip -c \
  | cut -f 2-4 \
  | get-regular-chromosomes \
  | sort-compress-index-bed --name ${gaps_path}/gaps.sorted



