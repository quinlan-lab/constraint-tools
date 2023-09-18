source set-environment-variables.sh 

# https://stackoverflow.com/a/43476575/6674256
# need to export PYTHONPATH since it is not already in the environment: 
# `printenv | grep -w PYTHONPATH` returns zero output
export PYTHONPATH="${CONSTRAINT_TOOLS}/utilities"

write-vista-elements () { 
  local label=$1 # "positive" or "negative"

  # https://enhancer.lbl.gov/
  url="https://enhancer.lbl.gov/cgi-bin/imagedb3.pl?search.form=no;search.result=yes;show=1;order=hs;page_size=100;action=search;form=search;page=1;search.sequence=1"

  output=${CONSTRAINT_TOOLS_DATA}/vista-elements/vista-elements.${label}.hg19.tsv

  curl $url 2> /dev/null \
    | python ${CONSTRAINT_TOOLS}/download-process-data/vista-elements/convert-to-bed.py ${label} 2> /dev/null \
    | sort --version-sort -k1,1 -k2,2n -k3,3n \
    > ${output}

  info "Wrote:" ${output}
}

write-vista-elements "positive"
write-vista-elements "negative"
