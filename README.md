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

Assuming one has access to the protected environment on the CHPC at University of Utah: 

```
bash tests/germline-model/train.sh $PWD 
```

Once training is complete, do: 
```
bash tests/germline-model/dashboard.sh $PWD
```

Follow the instructions at the command line to view a web app that visualizes observed SNV and singleton counts, and those expected under a null model of sequence-dependent mutation (see `define-model` folder), as a function of genomic coordinate.  

A plot of estimated SNV probabilities of the neutral model can be found here: https://github.com/quinlan-lab/constraint-tools/blob/main/tests/plot_mutation_probabilities.ipynb
 
## Usage

```
./constraint-tools [SUB_COMMAND] REQUIRED_ARGUMENTS
```

Valid values for `SUB_COMMAND` are: 

```
train-germline-model 
      estimate kmer-dependent SNV probabilities and singleton-count probabilities (see the model defined in the "define-model" folder)
dashboard-germline-model
      start a web app that visualizes observed and expected SNV and singleton counts as a function of genomic coordinate
predict-germline-constraint
      call genomic regions predicted to be under negative selection in the germline [not yet implemented]
```

Required arguments for `train-germline-model` are:

```
--genome STR
      Path to a reference fasta. 
      A "samtools faidx" index is expected to be present at the same path. 
--mutations STR 
      Path to a set of mutations specified as tab-separated values with column headings: "chromosome start end variant_type REF ALT number_ALT number_ALT_chromosomes number_chromosomes SYMBOL Gene Amino_acids CANONICAL Consequence Feature_type Feature miscellaneous". A "tabix" index is expected to be present at the same path.
--kmer-size INT
      Size of kmer of model to be trained. 
--model STR 
      Path to a directory to store trained model in. 
```

By default the `train-germline-model` subcommand uses a pre-computed set of putatively neutral regions from the GRCH38 reference located in the `/dist` folder, and a reasonable value of the window-size to compute singleton counts within. Optionally, the user may change either of these defaults by specifying the `--neutral-regions` and `--window-size` arguments: 

```
--neutral-regions STR
      Bed-format file containing a list of genomic intervals on which the model is to be trained.
--window-size INT
      size of the intervals used to compute the null distribution of singleton count.
```

Running `train-germline-model` produces a specification of the sequence-dependent and allele-frequency-aware neutral model in json format, viewable using, e.g., 
```
${CONSTRAINT_TOOLS}/bin/jq . ${model_folder}/<json file> 
```

Required arguments for `dashboard-germline-model` are:

```
--port INT 
      The port to serve the web-app on
```

By default the `dashboard-germline-model` subcommand uses a pre-computed model. 
Optionally, the user may change this by specifying the `--model` argument: 

```
--model STR
      Path to a neutral model produced by the train-germline-model sub-command (in json format). This model is used to compute the expected SNV and singleton counts in the visualization. 
```

## Input Data

Assuming one has access to the protected environment on the CHPC at University of Utah, 
then data files can be found at: 

```
/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools
```

## Production model

In the `/dist` directory, we distribute a model 
that was trained on a genome-wide set of putatively neutral regions
(also located in the `/dist` directory).

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
