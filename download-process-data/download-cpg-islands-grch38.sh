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

info "mysql path is:" "$(which mysql)"

cpg_islands_path="${CONSTRAINT_TOOLS_DATA}/cpg-islands/grch38"
mkdir --parents ${cpg_islands_path}

database="hg38"
table="cpgIslandExt"

# https://genome.ucsc.edu/goldenPath/help/mysql.html
# https://genome.ucsc.edu/cgi-bin/hgTables
mysql \
    --user=genome \
    --host=genome-mysql.soe.ucsc.edu \
    --port=3306 \
    --batch \
    --no-auto-rehash \
    --skip-column-names \
    --execute="
      SELECT 
      chrom, chromStart, chromEnd, name, length, cpgNum, gcNum, perCpg, perGc, obsExp 
      FROM ${table};
    " \
    ${database} \
  | get-regular-chromosomes \
  | sort-compress-index-bed --name ${cpg_islands_path}/cpg-islands.sorted
