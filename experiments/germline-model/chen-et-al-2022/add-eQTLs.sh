set -o errexit
set -o pipefail
# set -o noclobber
# set -o xtrace
set -o nounset

source set-environment-variables.sh 

PATH="${CONSTRAINT_TOOLS}/experiments/germline-model/chen-et-al-2022:$PATH" 

# https://stackoverflow.com/a/43476575/6674256
# need to export PYTHONPATH since it is not already in the environment: 
# `printenv | grep -w PYTHONPATH` returns zero output
export PYTHONPATH="${CONSTRAINT_TOOLS}/utilities"

chen_mchale_windows="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/chen-mchale.bed"

# eQTLs 
# download-process-data/eQTL/download.sh
eQTLs="/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/eQTL/figS40c_input_data.txt"

get-chen-windows () { 
  cat ${chen_mchale_windows} | tail -n +2
}

get-unique-eQTLs () {
  get-eQTLs.py ${eQTLs} \
    | sort --version-sort -k1,1 -k2,2 \
    | uniq 
}


# get-chen-windows | head 
# get-unique-eQTLs 
# exit 1 

header-line () {
  echo -e "chromosome\tchen_start\tchen_end\tchen_zscore\tmchale_start\tmchale_end\tmchale_position\tmchale_N_bar\tmchale_N_observed\tmchale_K_bar\tmchale_K_observed\tmchale_M\tchen_mchale_overlap_bps\teQTL_chromosome\teQTL_start\teQTL_end\twatershed_posterior"
}

add-eQTLs () {
  bedtools intersect \
      -a <(get-chen-windows) \
      -b <(get-unique-eQTLs) \
      -wa \
      -wb 
}

chen_mchale_windows_with_eQTLs="${CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/chen-mchale.eQTL.bed"
(
  header-line 
  add-eQTLs
) > ${chen_mchale_windows_with_eQTLs}  
# | head | column -t -s $'\t'

info "Wrote chen-mchale windows, with eQTLs, to" ${chen_mchale_windows_with_eQTLs}  

