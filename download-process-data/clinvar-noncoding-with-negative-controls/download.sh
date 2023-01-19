set -o errexit
set -o pipefail
set -o xtrace
set -o nounset

source download-process-data/set-environment-variables.sh 

# https://github.com/melobio/LOGO/blob/master/05_LOGO_Variant_Prioritization/2.%20data/4.%20Clinvar%20+%201000G_177Pos_177Neg/Clinvar_nc_snv_pathogenic_177Pos_177Neg.vcf 
clinvar_url='https://raw.githubusercontent.com/melobio/LOGO/master/05_LOGO_Variant_Prioritization/2.%20data/4.%20Clinvar%20%2B%201000G_177Pos_177Neg/Clinvar_nc_snv_pathogenic_177Pos_177Neg.vcf'
data_directory="/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/clinvar-noncoding-with-negative-controls"

curl ${clinvar_url} > $data_directory/Clinvar_nc_snv_pathogenic_177Pos_177Neg.hg19.vcf

# "Converting single nucleotide variants between genome builds: from cautionary tale to solution" : 
# https://academic.oup.com/bib/article/22/5/bbab069/6210068 
liftover_unstable_positions_url="https://raw.githubusercontent.com/cathaloruaidh/genomeBuildConversion/master/CUP_FILES/FASTA_BED.ALL_GRCh37.novel_CUPs.bed"

curl $liftover_unstable_positions_url > $data_directory/FASTA_BED.ALL_GRCh37.novel_CUPs.bed
