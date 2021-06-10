# constraint-tools

Mutation Annotation Format (MAF) specification: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/


## Gotchas

Based upon: https://www.biostars.org/p/69222/

1. Ensure that different tumor sample barcodes were not used for the same sample by reducing the tumor IDs to the form "TCGA-XX-XXXX-XX". 
2. Our model assumes that the probability of generating a mutation at a given site in a given tumor is independent of whether a mutation has been seen at the same site in another tumor. That assumption is violated when tumors share a common ancestor, as is the case for a primary tumor and a metastasis. Therefore one would probably want to ensure that no more than one tumor per patient is included in the maf.

 
