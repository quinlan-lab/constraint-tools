import sys
import numpy as np
import gzip
import subprocess
from random import sample

## A script that generates simulated bed files
## Takes a source file and shuffles those events into the genome
## Uses bedtools shuffle 
## Resulting events should match the originals, but with new locations 
## To be used as input for [simulation_pipeline].sh

## Hard coding a genome file so this isn't needed anymore
## Get chromosome sizes from genome file
#genome = open('/uufs/chpc.utah.edu/common/HIPAA/u0055382/genome_ref/GRCh38.frankenstein.genome', 'r')
#chrom_sizes = {}
#for line in genome:
#    chrom,size = line.rstrip('\n').split('\t')
#    if chrom not in chrom_sizes:
#        chrom_sizes[chrom] = size

## Designate an exclude and genome file for bedtools shuffle
exclude_file = '/scratch/ucgd/lustre-work/quinlan/u0055382/SV_constraint/exclude.bed'
genome_file = '/scratch/ucgd/lustre-work/quinlan/u0055382/SV_constraint/GRCh38.autosomes_only.genome'

## Load source file to base simulated data off of
## Save all original columns for later
## Save a list of chromosomes in the source file
## Generate tmp input file for bedtools shuffle
source_file = gzip.open(sys.argv[1], 'rt', encoding='utf-8')
saved_info = {}
chroms = []

## Need a unique identifier for tmp files
tmp_id = sys.argv[2]
tmp_input_name = 'tmp_input.' + str(tmp_id) + '.bed'
tmp_input = open(tmp_input_name, 'w')
for line in source_file:
    fields = line.rstrip('\n').split('\t')
    chrom = fields[0]

    ## Skipping events on X and Y chromosomes alternate ones with '_'
    ## Otherwise bedtools shuffle inflates the number of events on these chromosomes
    if chrom == 'X' or chrom == 'Y' or '_' in chrom:
        continue

    svlen = fields[3]
    sv_id = fields[6]

    if sv_id not in saved_info:
        saved_info[sv_id] = fields

    if chrom not in chroms:
        chroms.append(chrom)

    tmp_input.write('\t'.join([str(chrom),str(0),str(svlen),str(sv_id)]))
    tmp_input.write('\n')

tmp_input.close()

## Hard coded a genome file so this isn't needed anymore
## Generate a tmp genome file for bedtools shuffle
#tmp_genome = open('tmp.genome', 'w')
#for chrom in chroms:
#    tmp_genome.write('\t'.join([str(chrom),str(chrom_sizes[chrom])]))
#    tmp_genome.write('\n')

#tmp_genome.close()

## Run bedtools shuffle to create new genomic locations for the events
## An exclude file prevents events in certain location
## Creates tmp_shuffled.bed
tmp_shuffled_name = 'tmp_shuffled.' + str(tmp_id) + '.bed'
tmp_shuffled = open(tmp_shuffled_name, 'w')
subprocess.run(["bedtools",\
                "shuffle",\
                "-i",tmp_input_name,\
                "-g",genome_file, \
                "-excl",exclude_file], \
                stdout=tmp_shuffled)

tmp_shuffled.close()

## Read in tmp_shuffled.bed and then add extra columns
## Print out the results
shuffled = open(tmp_shuffled_name, 'r')
for line in shuffled:
    chrom,start,end,sv_id = line.rstrip('\n').split('\t')
    svlen = saved_info[sv_id][3]
    svtype = saved_info[sv_id][4]
    source = saved_info[sv_id][5]
    af = saved_info[sv_id][7]
    hom_ref = saved_info[sv_id][8]
    het = saved_info[sv_id][9]
    hom_alt = saved_info[sv_id][10]
    out = [chrom,start,end,svlen,svtype,source,sv_id,af,hom_ref,het,hom_alt]
    print('\t'.join(list(map(str,out))))
