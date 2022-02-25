Download gnomad vcfs: 
```
cd ${CONSTRAINT_TOOLS} 
bash download-process-data/gnomad/v3/variants/download-gnomad-variants-main.sh
```

Convert vcfs to tsvs (the format expected by `constraint-tools`): 

```
sbatch download-process-data/gnomad/v3/variants/vcf-to-tsv-chromosomes.sh
sbatch download-process-data/gnomad/v3/variants/merge-chromosomes.sh 
```