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

# # http://uswest.ensembl.org/info/docs/api/core/core_tutorial.html
# mysql --user=anonymous --host=ensembldb.ensembl.org --execute "show databases;" \
#   | grep homo_sapiens_core_

genes_path="${CONSTRAINT_TOOLS_DATA}/genes/grch38"
mkdir --parents ${genes_path}

# https://uswest.ensembl.org/info/data/mysql.html
# https://stackoverflow.com/a/16592892/6674256
mysql \
    --user=anonymous \
    --host=ensembldb.ensembl.org \
    --port=3306 \
    --no-auto-rehash \
    --batch \
    --skip-column-names \
    homo_sapiens_core_98_38 \
  < ${CONSTRAINT_TOOLS}/download-process-data/canonical-exons-grch38.sql \
  | get-regular-chromosomes \
  | awk '{ print "chr"$0 }' \
  | sort-compress-index-bed --name ${genes_path}/canonical-exons.sorted

