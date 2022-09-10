source set-environment-variables.sh 

# https://stackoverflow.com/a/43476575/6674256
# need to export PYTHONPATH since it is not already in the environment: 
# `printenv | grep -w PYTHONPATH` returns zero output
export PYTHONPATH="${CONSTRAINT_TOOLS}/utilities"

# https://enhancer.lbl.gov/
url="https://enhancer.lbl.gov/cgi-bin/imagedb3.pl?search.form=no;search.result=yes;show=1;order=hs;page_size=100;action=search;form=search;page=1;search.sequence=1"

curl $url \
  | grep "^>Human" \
  | python ${CONSTRAINT_TOOLS}/download-process-data/vista-enhancers/convert-to-bed.py \
  | sort --version-sort -k1,1 -k2,2n -k3,3n \
  > ${CONSTRAINT_TOOLS_DATA}/vista-enhancers/vista-enhancers.hg19.tsv