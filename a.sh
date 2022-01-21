source download-process-data/set-environment-variables.sh

# https://stackoverflow.com/a/43476575/6674256
# need to export PYTHONPATH since it is not already in the environment:
# `printenv | grep -w PYTHONPATH` returns zero output
export PYTHONPATH="${CONSTRAINT_TOOLS}/utilities"

# no need to export PATH since it is already in the environment:
# `printenv | grep -w PATH` returns non-zero output
PATH="${CONSTRAINT_TOOLS}/bin:${CONSTRAINT_TOOLS}/download-process-data/gnomad/v3/variants:$PATH"

vcf-to-tsv-interval \
  --interval chr1:5974969-5975969 \
  --vcf /scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/gnomad/v3/variants/gnomad_v3_chr1.vcf.gz \
  --vep-keys /scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/gnomad/v3/variants/gnomad_v3_chr1.vep-keys.txt \
  --tmpdir /scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/gnomad/v3/variants/tmp.aWpwTX7hGU
