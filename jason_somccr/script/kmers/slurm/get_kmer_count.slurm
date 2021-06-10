#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --nodes=4
#SBATCH --ntasks=16
#SBATCH --account=quinlan
#SBATCH --partition=redwood-shared-freecycle
#SBATCH -o get_kmer_count.%j.log

## Define variables
dir=$HOME/git/somccr

## Load modules
module use $HOME/MyModules
module load miniconda3/latest

#-----------------------------------------------------------------------------#

## Parallelize function to download raw hg19 UCSC fasta files via threading
python $dir/script/kmers/python/get_chr_fasta_files.py


#-----------------------------------------------------------------------------#

## Get unmutated intervals from the consensus variant file 
python $dir/script/kmers/python/get_constrained_region.py

#-----------------------------------------------------------------------------#

## Parallelize function to obtain kmer counts across all CDS intervals and variants (SNVs)
python $dir/script/kmers/python/get_kmer_count.py 

#-----------------------------------------------------------------------------#















