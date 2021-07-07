# constraint-tools

## Installation

```
git clone https://github.com/quinlan-lab/constraint-tools
cd constraint-tools
bash install.sh 
bash build-vue-app.sh
```
Only installation on Linux x86_64 is currently supported.

## Quick Start 

Assuming one has access to the protected environment on the CHPC at University of Utah: 

```
bash tests/test.sh
```

This starts a web app. Visit `localhost:5000` in your web browser to view the web app. 

If the software is run on a remote machine, then you can view the web-app on a local machine by doing the following from the local machine (change username and host first):

```
ssh -N -L localhost:5000:localhost:5000 username@host
```

A plot of estimated mutation probabilities that are fed into the model can be found here: https://github.com/quinlan-lab/constraint-tools/blob/main/tests/plot_mutation_probabilities.ipynb
 
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
visualize
      start a web app that visualizes mutation counts as a function of genomic coordinate
predict
      call genomic regions predicted to be under negative selection [not yet implemented]
```

Required arguments for `train` are:

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

This produces a specification of the sequence-dependent mutation model in json format, viewable using, e.g., 
```
${root}/bin/jq . ${output}/<json file> 
```

Required arguments for `visualize` are:

```
--model STR
      Path to the model produced by the train sub-command (in json format)
--port INT 
      The port to serve the web-app on
```
      
## Input 

Assuming one has access to the protected environment on the CHPC at University of Utah, 
then sorted, block-compressed, and indexed vcf, maf, gtf and fasta files can be found at: 

```
/scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools/data
```

## Development 

Changes to `vue-app` necessitate a rebuilding of the vue app by running 

```
bash build-vue-app.sh 
```

## TODO 

Jason will manually select a contiguous, putatively neutral region of the genome to train the model on using
`constraint-tools train ...`. He will also find the genomic coordinates of a positive-control exon. 

Peter will implement a minimal version of `constraint-tools visualize ...`. 

## Mutation Annotation Format (MAF) 

1. specification: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/
2. minimal example: https://github.com/mskcc/vcf2maf/blob/main/data/minimalist_test_maf.tsv
3. tooling: 
    1. [maf-lib](https://github.com/NCI-GDC/maf-lib): comprehensive, but lacks so much documentation that it is effectively unusable
    2. [vcf2maf](https://github.com/mskcc/vcf2maf): can convert maf to vcf, but reliance on vep makes tool effectively unusable
4. gotchas: https://www.biostars.org/p/69222/


