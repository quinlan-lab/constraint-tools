```
cd ${CONSTRAINT_TOOLS} 
sbatch download-process-data/gnomad/v3/variants/vcf-to-tsv-chromosomes.sh
sbatch merge-chromosomes.sh 
```