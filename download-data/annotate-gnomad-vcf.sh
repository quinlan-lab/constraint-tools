#!/bin/sh
#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --nodes=2
#SBATCH --ntasks=16
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH -o logs/annotate-gnomad-vcf.out

source download-data/set-environment-variables.sh 
module load bcftools
module load vep

## Define directory to download files into 
gnomad_path="${CONSTRAINT_TOOLS}/data/gnomad/v3"
mkdir --parents ${gnomad_path}

## Define variables 
gnomad_variant="gnomad_v3_chr22.vcf.bgz"
gnomad_tbi="gnomad_v3_chr22.vcf.bgz.tbi"
gnomad_variant_filtered_sorted="gnomad_v3_chr22_filtered_sorted.vcf.bgz"
gnomad_variant_filtered_sorted_annotated="gnomad_v3_chr22_filtered_sorted.vep.gz"

info "Filtering gnomad vcf to remove indels..."
cd ${gnomad_path}
bcftools view gnomad_v3_chr22.vcf.bgz -H --types snps | sort -k1,1 -k2,2n | bgzip -c  > ${gnomad_variant_filtered_sorted}

info "Annotating gnomad vcf with VEP..."
vep -i ${gnomad_variant_filtered_sorted} -o ${gnomad_variant_filtered_sorted_annotated} --everything --compress_output bgzip







