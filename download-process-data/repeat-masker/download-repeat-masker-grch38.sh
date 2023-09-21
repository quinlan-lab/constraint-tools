set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset 

source download-process-data/set-environment-variables.sh 

# no need to export PATH since it is already in the environment: 
# `printenv | grep -w PATH` returns non-zero output 
PATH="${CONSTRAINT_TOOLS}/bin:$PATH" 
PATH="${CONSTRAINT_TOOLS}/utilities:$PATH" 
PATH="${CONSTRAINT_TOOLS}/download-process-data/repeat-masker:$PATH" 

info "mysql path is:" "$(which mysql)"

repeat_masker_path="${CONSTRAINT_TOOLS_DATA}/repeat-masker/grch38"
mkdir --parents ${repeat_masker_path}

database="hg38"
table="rmsk" # https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=1713033924_9X9h4NmBdDKuvdwCDL0ACgaOyMV7&db=hg38&g=rmsk
column_names="genoName, genoStart, genoEnd, strand, repName, repClass, repFamily" 

# schema: 
# https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=rep&hgta_track=rmsk&hgta_table=rmsk&hgta_doSchema=describe+table+schema

# https://genome.ucsc.edu/goldenPath/help/mysql.html
mysql \
    --user=genome \
    --host=genome-mysql.soe.ucsc.edu \
    --port=3306 \
    --batch \
    --no-auto-rehash \
    --skip-column-names \
    --execute="
      SELECT 
      ${column_names}
      FROM ${table};
    " \
    ${database} \
  | get-regular-chromosomes \
  | get-SINE-LINE-LTR-DNA \
  | sort-compress-index-bed --name ${repeat_masker_path}/repeat-masker.SINE-LINE-LTR-DNA.sorted

echo "${column_names}" > ${repeat_masker_path}/repeat-masker.SINE-LINE-LTR-DNA.column-names.csv

