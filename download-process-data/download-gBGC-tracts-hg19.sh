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

gBGC_path="${CONSTRAINT_TOOLS_DATA}/GC-biased-gene-conversion"
mkdir --parents ${gBGC_path}

database="hg19"
table="phastBiasTracts3" # https://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=compGeno&hgta_track=phastBias&hgta_table=phastBiasTracts3&hgta_doSchema=describe+table+schema
column_names="chrom, chromStart, chromEnd" 

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
  | sort --parallel=8 --buffer-size=75% --version-sort -k1,1 -k2,2 \
  > ${gBGC_path}/gBGC-tracts.hg19.sorted.bed

echo "${column_names}" > ${gBGC_path}/gBGC-tracts.hg19.column-names.csv

lift () {
  local filename_root=$1
  info "Lifting ${filename_root}.hg19.sorted.bed to hg38"
  bash \
    ${CONSTRAINT_TOOLS}/download-process-data/lift.sh \
    "${filename_root}.hg19.sorted.bed" \
    hg19 \
    hg38
  mv \
    ${filename_root}.hg19.sorted.bed.hg38 \
    ${filename_root}.hg38.bed
  mv \
    ${filename_root}.hg19.sorted.bed.hg38.unmapped \
    ${filename_root}.hg19.unmapped.bed
}

lift "${gBGC_path}/gBGC-tracts"

