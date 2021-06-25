# constraint-tools

## Installation

```
git clone https://github.com/quinlan-lab/constraint-tools
cd constraint-tools
bash install.sh 
```
Only installation on Linux x86_64 is currently supported.

## Quick Start 

Assuming one has access to the protected environment on the CHPC at University of Utah: 

```
bash tests/test.sh
```

Then view https://github.com/quinlan-lab/constraint-tools/blob/main/tests/plot_mutation_probabilities.ipynb
 
## Usage

Assuming that the path to the `constraint-tools` directory on your filesystem is `${root}`, usage is:

```
PATH="${root}:$PATH"
constraint-tools [SUB_COMMAND] REQUIRED_ARGUMENTS
```

Valid sub-commands are: 

```
train 
      estimate kmer-dependent mutation probabilities 
predict
      call genomic regions predicted to be under negative selection 
```

Required arguments are:

```
--genome STR
      Path to the reference fasta. 
      A "samtools faidx" index is expected to be present at the same path. 
--region STR 
      Samtools-style specification of a genomic interval (see examples in tests directory).
--mutations STR 
      Path to a set of mutations specified in Mutation Annotation Format.
      A "tabix" index is expected to be present at the same path.
--kmer-size INT
      Size of kmer to use in model. 
--output STR 
      Path to a directory to store results in. 
```

## Input 

Assuming one has access to the protected environment on the CHPC at University of Utah, 
then sorted, block-compressed, and indexed vcf, maf, gtf and fasta files can be found at: 

```
/scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools/data
```

## Output 

A specification of the sequence-dependent mutation model in json format, viewable using, e.g., 
```
${root}/bin/jq . <json file> 
```

## TODO 

Jason will manually select a contiguous, putatively neutral region of the genome to train the model on using
`constraint-tools train ...`. He will also find the genomic coordinates of a positive-control exon. 

Peter will implement a minimal version of `constraint-tools predict ...`. 

## Mutation Annotation Format (MAF) 

1. specification: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/
2. minimal example: https://github.com/mskcc/vcf2maf/blob/main/data/minimalist_test_maf.tsv
3. tooling: 
    1. [maf-lib](https://github.com/NCI-GDC/maf-lib): comprehensive, but lacks so much documentation that it is effectively unusable
    2. [vcf2maf](https://github.com/mskcc/vcf2maf): can convert maf to vcf, but reliance on vep makes tool effectively unusable
4. gotchas: https://www.biostars.org/p/69222/


