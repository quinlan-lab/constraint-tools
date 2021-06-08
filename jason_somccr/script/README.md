# **Description**
This file describes the workflow to identify regions of constraint in tumor genomes. Specifically, I seek to identify intervals of coding sequences in tumor genomes signficantly devoid of protein-altering mutations. This analysis will be done at a pan-cancer and individual cancer type analysis. 

We will use consensus variant files from **TCGA (whole exome sequencing)** and **ICGC (whole genome sequencing)** working groups. See below for more information on the working groups and their manuscripts: 

1.   [Scalable Open Science Approach for Mutation Calling 
of Tumor Exomes Using Multiple Genomic Pipelines](https://pubmed.ncbi.nlm.nih.gov/29596782/)
- Use publically available MAF file: [mc3.v0.2.8.PUBLIC.maf.gz](https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc)
- Coverage information:  
- Clinical information: Download a tsv file containing TCGA clinical data from the [NCI GDC Data Portal](https://portal.gdc.cancer.gov/repository?searchTableTab=cases) by clicking on the **"Clinical"** button. 

```
## Decompress TCGA clinical information file:
cd ~git/somccr/data/clinical
tar -zxvf ~/path/to/downloaded/file > mc3_mapping.tsv
```
    
2.   [Pan-cancer analysis of whole genomes](https://www.nature.com/articles/s41586-020-1969-6) 
* Use publically available MAF file: [final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz](https://dcc.icgc.org/api/v1/download?fn=/PCAWG/consensus_snv_indel/final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz)
* Use publically available wig files for coverage information: [coverage_wig_files.tar](https://dcc.icgc.org/api/v1/download?fn=/PCAWG/consensus_snv_indel/wig_files/coverage_wig_files.tar)
* Clinical information: Click the **"Download Donor Data"** button at [ICGC data portal](https://dcc.icgc.org/search?filters=%7B%22donor%22:%7B%22id%22:%7B%22is%22:%5B%22ES:12b6fcab-467d-4649-8177-0aa41f44d77c%22%5D%7D%7D%7D). Use the **submitted_donor_id** and **project_code** columns.  
    
```
## Decompress ICGC clinical information file:
cd ~git/somccr/data/clinical
gunzip ~/Path/to/sample.tsv.gz > pcawg_mapping.tsv
```

# Import Directories: 

* `/scratch/ucgd/lustre-work/quinlan/u1240855/somccr/data` --> stores large files
    * output
    * reference
        * cds
        * chr
        * segdup_selfchain
* `/uufs/chpc.utah.edu/common/HIPAA/u1240855/git/somccr`
    * data --> stores ICGC and TCGA clinical information
    * presentations
    * scripts


# **Step 1: Format, Filter, and Combine TCGA and ICGC Variants**
## ***Step 1a: Obtain variant and coverage files from working groups***
### Step 1a.1: Get variant information from TCGA and ICGC 
### Step 1a.2: Map clinical information (i.e. cancer type) to sample ids from both working group datasets

### **Generate consensus cancer type mapping file** 
* **GOAL:** Generate a consensus file to map ICGC's 'project_code' to MC3's 'project_id' to the same disease 
* Use [TCGA study abbreviations](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations) and [ICGC study abbreviations](https://docs.icgc.org/submission/projects/) to obtain consensus TCGA and ICGC cancer type name. 
* The consensus cancer type name will be indicated in the **'Simplified Name'** column...

| TCGA Study Code | ICGC Project Code | ICGC Cancer Type | TCGA Cancer Type |Simplified Name | TCGA? | ICGC? |
| --------------- | ----------------- | ---------------- | ---------------- | -------------- | ----- | ----- |
| TCGA-BRCA | BRCA-CN | Breast Triple Negative Cancer |  Breast invasive carcinoma |Breast invasive carcinoma | 1 | 1 |
| TCGA-BRCA | BRCA-EU | Breast ER+ and HER2 | Breast invasive carcinoma | Breast invasive carcinoma | 1 | 1 |
| TCGA-BRCA | BRCA-FR | Breast Cancer | Breast invasive carcinoma | Breast invasive carcinoma | 1 | 1 |
| TCGA-BRCA | BRCA-KR | Breast Cancer | Breast invasive carcinoma | Breast invasive carcinoma | 1 | 1 |
| TCGA-BRCA | BRCA-UK | Breast Triple Negative/Lobular Cancer | Breast invasive carcinoma | Breast invasive carcinoma | 1 | 1 |

* **Some cancer type names are only found in ICGC --> see below for a few examples:**

| TCGA Study Code | ICGC Project Code | ICGC Cancer Type | TCGA Cancer Type |Simplified Name | TCGA? | ICGC? |
| --------------- | ----------------- | ---------------- | ---------------- | -------------- | ----- | ----- |
| NA | ALL-US | Acute Lymphoblastic Lymphoma | NA | Acute Lymphoblastic Lymphoma | 0 | 1 |
| NA | LIAD-US | Benign Liver Tumor | NA | Benign Liver Tumor | 0 | 1 |
| NA | CCSK-US | Clear Cell Sarcomas of the Kidney | NA | Clear Cell Sarcomas of the Kidney | 0 | 1 |
| NA | PEME-CA | Pediatric Medulloblastoma | NA | Pediatric Medulloblastoma | 0 | 1 |
| NA | RT-US | Rhabdoid Tumor | NA | Rhabdoid Tumor | 0 | 1 |
| NA | SKCA-BR | Skin Adenocarcinoma | NA | Skin Adenocarcinoma | 0 | 1 |
| NA | LMS-FR | Soft tissue cancer | NA | Soft tissue cancer | 0 | 1 |
| NA | WT-US | Wilms Tumor | NA | Wilms Tumor | 0 | 1 |

* **Some cancer type names are only found in TCGA --> see below for a few examples:**

| TCGA Study Code | ICGC Project Code | ICGC Cancer Type | TCGA Cancer Type |Simplified Name | TCGA? | ICGC? |
| --------------- | ----------------- | ---------------- | ---------------- | -------------- | ----- | ----- |
| TCGA-UVM | NA | NA | Uveal Melanoma | Uveal Melanoma | 1 | 0 |
| TCGA-MESO | NA | NA | Mesothelioma | Mesothelioma | 1 | 0 |
| TCGA-THYM | NA | NA | Thymoma | Thymoma | 1 | 0 |

### Step 1a.3 Filter variants based on coverage parameters

## ***Step 1b: Remove variants that overlap with regions of segmental duplications and/or self-chains***

### Step 1b.1: Concatenate variants from MC3 and PCAWG MAF files

### Step 1b.2: Select for specific variant classes (exonic and/or non-intronic)

### Step 1b.3: Use bedtools to remove variants that overlap with **regions of segdups and self-chains** (see `/scratch/ucgd/lustre-work/quinlan/u1240855/somccr/data/reference/segdup_selfchain/get_segdup_selfchain.ipynb`)

### Step 1b.4: Use bedtools to select for variants within CDS regions of genes
#### Get CDS gene coordiantes from UCSC (hg19)

* **Store the file in `~/git/somccr/data/reference/segdup_selfchain`**

Table Browser             |  Output Query
:-------------------------:|:-------------------------:
![Table Browser](https://gitlab.chpc.utah.edu/u1240855/somccr/-/blob/master/presentations/img/table_browser_query.png)  |  ![Output Query](https://gitlab.chpc.utah.edu/u1240855/somccr/-/blob/master/presentations/img/get_output_query.png)

#### Use [bedtools merge](https://bedtools.readthedocs.io/en/latest/content/tools/merge.html) to combine sorted gene CDS intervals
```
cd  ~/git/somccr/data/reference/cds
sort -k1,1 -k2,2n cds_ucsc_gene_coord.bed | bedtools merge > cds_ucsc_gene_coord.sorted.merged.bed
cat cds_ucsc_gene_coord.sorted.merged.bed | sed 's/chr//g' > cds_ucsc_gene_coord.sorted.merged.noprefix.bed
cat cds_ucsc_gene_coord.sorted.merged.noprefix.bed | awk 'BEGIN {FS=OFS="\t"}; $1=="X" {$1=23}; {print $0}' | awk 'BEGIN {FS=OFS="\t"}; $1=="Y" {$1=24}; {print $0}' > cds_ucsc_gene_coord.sorted.merged.noprefix.renameXY.bed
cat cds_ucsc_gene_coord.sorted.merged.noprefix.renameXY.bed | awk '/^[0-9]*\t/ {print $0}' > cds_ucsc_gene_coord.sorted.merged.noprefix.renameXY.onlyChr.bed
```

### Step 1b.5: Use biomaRt to filter for genes with known CDS lengths (using hg19) 

# **Step 2: Calculate k-mer Frequencies**

* **See [k-mer counting tutorial](https://bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction/) for more information**

Reference k-mer | Mutated k-mer | Final Label | Frequency
:-: | :-: | :-: | :-:
ATC | AAC | ATC > AAC | ##%
ATC | ACC | ATC > ACC | ##%
ATC | AGG | ATC > AGC | ##%
## Step 2a: Get k-mer counts in the hg19 reference genome from all CDS regions
## Step 2b: Get k-mer counts for mutated k-mers
## Steb 2c: Get alternate k-mer frequencies

