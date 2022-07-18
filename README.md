# constraint-tools

## Installation

```
git clone https://github.com/quinlan-lab/constraint-tools
cd constraint-tools
bash install.sh 
bash build-vue-app.sh
```

This creates a conda environment that should be activated with, e.g.,:
```
conda activate constraint-tools
```

Only installation on Linux x86_64 is currently supported. 
Tested in the Protected Environment computer cluster of the Center for High Performance Computing (CHPC) at University of Utah. 

## Quick Start 

Assuming one has access to the protected environment on the CHPC at University of Utah,
one can do a prototype training using: 

```
bash tests/germline-model/train-germline-model-fast.sh
```

One can render a prototype browser using: 
```
bash tests/germline-model/browse.sh
```

Follow the instructions at the command line to view a web app that visualizes observed SNV and singleton counts, and those expected under a null model of sequence-dependent mutation (see `define-model` folder), as a function of genomic coordinate.  
 
## Usage

```
./constraint-tools [SUB_COMMAND] REQUIRED_ARGUMENTS
```

Valid values for `SUB_COMMAND` are: 

```
train-germline-model 
      Estimate kmer-dependent SNV probabilities and singleton-count probabilities 
      (see the model defined in the "define-model" folder)
browse-germline-model
      Start a web app that visualizes observed and expected 
      SNV and singleton counts as a function of genomic coordinate
predict-germline-model
      Tile trustworthy part of noncoding genome and compute z-scores for each tile. 
```

Required arguments for `train-germline-model` are:

```
--genome STR
      Path to a reference fasta. 
      A "samtools faidx" index is expected to be present at the same path. 
--build STR 
      Human reference genome build. 
      Allowed values for STR are "hg19" and "hg38".
--mutations STR 
      Path to a set of mutations specified as tab-separated values with column headings: 
      "chromosome start end variant_type REF ALT number_ALT 
      number_ALT_chromosomes number_chromosomes SYMBOL Gene Amino_acids 
      CANONICAL Consequence Feature_type Feature miscellaneous". 
      A "tabix" index is expected to be present at the same path.
--number-chromosomes-min INT
      Only consider SNVs at which the nucleotide identity (allele) is known 
      in >=INT chromosomes in the cohort.
--kmer-size INT
      Size of kmer in model to be trained. 
--model STR 
      JSON file to store the trained model in. 
--work STR 
      Path to a directory to store intermediate work and logs.
--progress-bars STR 
      Allowed values are "disk" or "stdout", 
      indicating whether to store logs containing "progress bars" 
      to disk or stdout, respetively.
```

By default the `train-germline-model` subcommand uses a pre-computed training set of trustworthy noncoding regions from the GRCH38 reference located in the `/dist` folder, and a reasonable value of the size of the window within which to count singletons. Optionally, the user may change either of these defaults by specifying the `--trustworthy-noncoding-regions` and `--window-size` arguments: 

```
--trustworthy-noncoding-regions STR
      Bed-format file containing a list of genomic intervals on which the model is to be trained.
--window-size INT
      Size of the intervals used to compute the null distribution of singleton count. 
      This is also the size of the window in "test" regions.
```

Other optional arguments are: 

```
--number-of-jobs INT 
      Number of slurm jobs to use during training. 
--max-trustworthy-noncoding-region-length INT 
      Trustworthy noncoding regions longer than this number are filtered out 
      from the set of regions that are ultimately used to train the model. 
```

Running `train-germline-model` produces a specification of the sequence-dependent and allele-frequency-aware null 
model in json format, viewable using, e.g., 
```
${CONSTRAINT_TOOLS}/bin/jq . <model> 
```

Required arguments for `browse-germline-model` are:

```
--port INT 
      The port to serve the web-app on
```

By default the `browse-germline-model` and `predict-germline-model` subcommands use a pre-computed model
corresponding to a pre-computed training set of trustworthy noncoding regions. 
Optionally, the user may change this by specifying the `--model` argument: 

```
--model STR
      Path to a null model produced by the train-germline-model sub-command 
      (in json format). 
      This model is used to compute the expected SNV and singleton counts. 

```

By default the `browse-germline-model` and `predict-germline-model` subcommands use a pre-computed set of trustworthy noncoding regions.
(A random subset of these regions comprise the default training set of trustworthy noncoding regions 
used to train the default model.)
Optionally, the user may change the set of trustworthy noncoding regions displayed 
by specifying the `--trustworthy-noncoding-regions` argument: 

```
--trustworthy-noncoding-regions STR
      Path to a set of trustworthy noncoding regions. 
      Ideally, of these, the trustworthy noncoding regions the model was trained upon should be a random subset. 
```

## Input Data

Assuming one has access to the protected environment on the CHPC at University of Utah, 
then data files can be found at: 

```
/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools
```

## Production models

In the `/dist` directory, we distribute models 
that were trained on a genome-wide training set of trustworthy noncoding regions
(also located in the `/dist` directory).
The bash script that was used to generate these models is:

```
train-germline-models-production-window-sizes.sh
```

## Sanity checks

See `experiments/germline-model/sanity-check-*.ipynb`. 

## Development note

Changes to the `vue-app` directory necessitate rebuilding the vue app by running 

```
bash build-vue-app.sh 
```

TODO: automate this using a git hook 

## Resources 
### Mutation Annotation Format (MAF) 

1. specification: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/
2. minimal example: https://github.com/mskcc/vcf2maf/blob/main/data/minimalist_test_maf.tsv
3. tooling: 
    1. [maf-lib](https://github.com/NCI-GDC/maf-lib): comprehensive, but lacks so much documentation that it is effectively unusable
    2. [vcf2maf](https://github.com/mskcc/vcf2maf): can convert maf to vcf, but reliance on vep makes tool effectively unusable
4. gotchas: https://www.biostars.org/p/69222/

### On k-mer counting 
https://bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction/
