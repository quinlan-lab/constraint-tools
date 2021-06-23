# constraint-tools

## Installation

git clone https://github.com/petermchale/trfermikit
cd trfermikit
bash install.sh 
conda activate trfermikit
Only installation on Linux x86_64 is currently supported.

## Data 

Sorted, block compressed, and indexed vcf, maf, gtf and fasta files can be found at: 

```
/scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools/data
```

## TODO 

1. get minimal CLI working, using a test set (one of Jason's positive control genes and a manually selected, small neutral region)
2. remove variants that have low read coverage, lie in low-complexity sequence, etc
3. stratify variants by cancer type 

## Mutation Annotation Format (MAF) 

1. specification: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/
2. minimal example: https://github.com/mskcc/vcf2maf/blob/main/data/minimalist_test_maf.tsv
3. tooling: 
    1. [maf-lib](https://github.com/NCI-GDC/maf-lib): comprehensive, but lacks so much documentation that it is effectively unusable
    2. [vcf2maf](https://github.com/mskcc/vcf2maf): can convert maf to vcf, but reliance on vep makes tool effectively unusable

## Data Gotchas

Based upon: https://www.biostars.org/p/69222/

1. Ensure that different tumor sample barcodes were not used for the same sample by reducing the tumor IDs to the form "TCGA-XX-XXXX-XX". 
2. Our model assumes that the probability of generating a mutation at a given site in a given tumor is independent of whether a mutation has been seen at the same site in another tumor. That assumption is violated when tumors share a common ancestor, as is the case for a primary tumor and a metastasis. Therefore one would probably want to ensure that no more than one tumor per patient is included in the maf.

