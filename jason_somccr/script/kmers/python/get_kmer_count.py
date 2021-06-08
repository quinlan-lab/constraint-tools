#!/usr/bin/env python3

"""
Obtain kmer counts from the reference genome
"""

#%%
from pyfaidx import Fasta
import numpy as np
import pandas as pd
import concurrent.futures
import collections 
from collections import defaultdict
import functools
import operator
import math
from scipy.stats import poisson

#%%
"""
Define function to count kmers from a given sequence
"""
def count_kmers(seq, kmer_size): 
    
    print('getting kmer counts in the interval')

    
    ## Initialize kmer count dictionary
    kCount = defaultdict(int)
    
    ## Iterate through each kmer available for a given sequence
    for i in range(0, len(seq) - kmer_size + 1):
        '''
        Get reference kmer counts: 
            See if kmer is already in the master kCount dictionary
        '''
        ## Get the i-th kmer
        kmer = str(seq[i:i + kmer_size])
        
        ## Add the kmer and kmer count to the dictionary
        if kmer in kCount:
            kCount[kmer] += 1
        else: 
            kCount[kmer] = 1
    
    return kCount

#%%
"""
Define overall function to count kmers within a given sequence
"""
def get_ref_kmer_counts(interval, kmer_size): 
    
    print(interval)
    
    ## Humanize variables
    chr = 0
    start = 1
    stop = 2
    
    ## Get the interval information
    chromosome = "chr" + str(int(interval[chr]))
    if chromosome == "chr23" :
        chromosome = "chrX"
    if chromosome == "chr24": 
        chromosome = "chrY"
    
    start = int(interval[start])
    stop = int(interval[stop])
    print(f'Defining the reference interval {chromosome}:{start}-{stop}')
    # return f'Defining the interval {chromosome}:{start}-{stop}'

    ## Get the fasta sequence of the interval
    print(f'getting fasta sequence for interval {chromosome}:{start}-{stop}')
    fasta_file = '/scratch/ucgd/lustre-work/quinlan/u1240855/somccr/data/reference/chr/' + chromosome + '.fa'
    fasta = Fasta(fasta_file, sequence_always_upper=True)
    seq = fasta.get_seq(chromosome, start, stop).seq
    
    '''
    Get kmer counts from the reference genome...
    '''
    ## Get the kmer counts for the interval's sequence
    kmer_count_dict = count_kmers(seq=seq, kmer_size=kmer_size) ## figure out how to allow this to be customizable
        
    return kmer_count_dict    


#%%
'''
Get reference_allele kmer from variant file --> assign it as dictionary key
Get alternate_allele from variant file --> assign it as value 
'''
def get_alt_kmer_counts(chromosome, start, stop, alt_allele, kmer_size):
    
    ## Get the interval information
    chromosome = "chr" + str(int(chromosome))
    if chromosome == "chr23" :
        chromosome = "chrX"
    if chromosome == "chr24": 
        chromosome = "chrY"
        
    start = int(start)
    stop = int(stop)
    alt_allele = str(alt_allele)
    
    ## Define the kmer interval
    start = start - kmer_size//2 + 1
    stop = stop + kmer_size//2
    print(f'Defining the alternate allele interval {chromosome}:{start}-{stop}...')

    ## Get the fasta sequence of the interval
    print(f'getting nucleotide sequence for interval {chromosome}:{start}-{stop}')
    fasta_file = '/scratch/ucgd/lustre-work/quinlan/u1240855/somccr/data/reference/chr/' + chromosome + '.fa'
    fasta = Fasta(fasta_file, sequence_always_upper=True)
    seq = fasta.get_seq(chromosome, start, stop).seq
        
    ## Get the reference kmer 
    ref_allele_kmer = seq
    
    ## Initialize dictionary with ref_allele_kmer as the key
    kmer_count_dict = {}
    kmer_count_dict[ref_allele_kmer] = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    
    ## Store alternate allele kmer
    kmer_count_dict[ref_allele_kmer][alt_allele] += 1
    
    return kmer_count_dict

#%%
"""
Pseudocode for reference kmer: 
1) Read in the fasta file
2) Iterate through each CDS bed interval via multiprocessing
3) Get the fasta sequence for each interval
4) Count kmers within each interval
5) Add all kmer counts across all intervals
    
Pseudocode for alternate kmer: 
1) 
"""

#%%
'''
Read in nonsilent mutations from the consensus variant file to get alternate kmer counts...
'''
## Read in the cds bed coordinates as a pandas dataframe
variant_file = '/scratch/ucgd/lustre-work/quinlan/u1240855/somccr/data/output/tcga_icgc.sorted.filtered.snv.CDS.knownCDSlength.no_chr_prefix.fix_header.renameXY.bed'
variant_df = pd.read_csv(variant_file, sep="\t", usecols=('chromosome', 'start', 'stop', 'variant_class', 'alt'))

## Select for non-silent mutations
nonsilent_df = variant_df.loc[(variant_df['variant_class'] != 'Silent') & (variant_df['alt'].isin(['A', 'C', 'G', 'T']))]

'''
Read in CDS regions to get reference kmer counts...
'''
bed_file = '/scratch/ucgd/lustre-work/quinlan/u1240855/somccr/data/reference/cds/cds_ucsc_gene_coord.sorted.merged.noprefix.renameXY.onlyChr.bed'
bed_regions = np.loadtxt(bed_file, delimiter="\t")

#%%
'''
Use multiprocessing to get reference kmer counts in regions of interest
'''
with concurrent.futures.ProcessPoolExecutor() as executor:
    ## Multiprocessing to get reference kmer counts
    kmer_size = [3] * int(bed_regions.shape[0])
    ref_results = executor.map(get_ref_kmer_counts, bed_regions, kmer_size)
    
    ## Combine kmer counts across all multiprocessing dictionaries
    my_list = []
    for result in ref_results: 
        my_list.append(result)
    
    ## Add kmer counts with the same kmer values
    combined_kmers = functools.reduce(operator.add, map(collections.Counter, my_list))
    
    ## Sort the reference kmers
    ref_kmers = collections.OrderedDict(sorted(combined_kmers.items()))
    
    print(ref_kmers)

#%%    
'''
Use multiprocessing to get alternate kmer counts from variant files
'''
with concurrent.futures.ProcessPoolExecutor() as executor:
    ## Multiprocessing to get alt kmer counts
    #nonsilent_df = nonsilent_df[:5]
    kmer_size = [3] * nonsilent_df.shape[0]
    alt_results = executor.map(get_alt_kmer_counts, nonsilent_df['chromosome'], nonsilent_df['start'],
                                    nonsilent_df['stop'], nonsilent_df['alt'], kmer_size)
       
    ## Combine kmer counts across all multiprocessing dictionaries
    my_list = []
    for result in alt_results: 
        my_list.append(result)
     
## Iterate through each variant kmer and combine the data
alt_kmer_dict = {}

for var_kmer in range(len(my_list)):
    
    ## Get the reference kmer
    ref_kmer = list(my_list[var_kmer].keys())[0]
    
    ## Get the alternate allele kmers
    alt_kmer = list(my_list[var_kmer].values())[0]
    
    ## Add the kmer and kmer count to the dictionary
    if ref_kmer in alt_kmer_dict: 
        alt_kmer_dict[ref_kmer] = functools.reduce(operator.add, map(collections.Counter, [alt_kmer_dict[ref_kmer], alt_kmer]))
    else: 
        alt_kmer_dict[ref_kmer] = alt_kmer
        

print(alt_kmer_dict)

#%%
'''
Function to process constrained intervals from the variant files: 
    1) Get sequnece of the constrained interval 
    2) Get the reference kmers present in the interval
    3) Get the frequency that each reference kmer was mutated to another nucleotide
    4) Sum up frequnecies 
    5) Get the length of the interval
    6) Calculate expected number of mutations
    7) Run Poisson probability for seeing 0 mutations with the expectation from #6
'''

def process_constrained_intervals(interval, ref_kmers, alt_kmer_dict, kmer_size):
    
    ## Humanize variables
    chr = 0
    start = 1
    stop = 2
    
    ## Get the interval information
    ccr_chromosome = "chr" + str(int(interval[chr]))
    if ccr_chromosome == "chr23" :
        ccr_chromosome = "chrX"
    if ccr_chromosome == "chr24": 
        ccr_chromosome = "chrY"
    
    ccr_start = int(interval[start])
    ccr_stop = int(interval[stop])
    
    ## Get the fasta sequence of the interval
    print(f'getting nucleotide sequence for interval {ccr_chromosome}:{ccr_start}-{ccr_stop}')
    fasta_file = '/scratch/ucgd/lustre-work/quinlan/u1240855/somccr/data/reference/chr/' + ccr_chromosome + '.fa'
    fasta = Fasta(fasta_file, sequence_always_upper=True)
    seq = fasta.get_seq(ccr_chromosome, ccr_start, ccr_stop).seq
    
    ## Initialize list to store frequencies of alterante kmer alleles
    alt_freq_list = []
    
    ## Iterate through each kmer in the constrained interval
    for i in range(0, len(seq) - kmer_size + 1):
        
        ## Get the i-th kmer
        kmer = str(seq[i:i + kmer_size])
        
        ## PLACEHOLDER --> REPLACE N WITH T
        kmer = kmer.replace("N", "T")
        
        ## Get the reference kmer count
        ref_kmer_count = ref_kmers[kmer]
        
        ## Get the alternate kmers
        alt_kmers = alt_kmer_dict[kmer]
        
        ## Get the sum of alternate allele kmers
        alt_kmer_count = sum(alt_kmers.values())
        
        ## Get the frequency of alternate kmers
        alt_frequency = alt_kmer_count / ref_kmer_count
        alt_freq_list.append(alt_frequency)
        
    ## Get the total mutant frequency over the interval of interest
    total_freq = sum(alt_freq_list)
    
    ## Get the interval length
    length = ccr_stop - ccr_start - 2
    
    ## Get the expected number of intervals
    expected_mut_num = math.ceil(total_freq * length)
    
    ## Get the Poisson probability
    prob = poisson.pmf(0, expected_mut_num)
    
    ## Store the interval information as a list
    interval_info = [ccr_chromosome, ccr_start, ccr_stop, expected_mut_num, prob]
    
    return interval_info

#%%
'''
Use multiprocessing to get seq of an interval of potential CCR
'''
## Get the constrained intervals file
constrained_intervals_file = '/scratch/ucgd/lustre-work/quinlan/u1240855/somccr/data/output/constrained_intervals/constrained_interval.bed'
constrained_intervals = np.loadtxt(constrained_intervals_file, delimiter="\t", skiprows=1)
#constrained_intervals = constrained_intervals[:5]

#%%
## Set up multiprocessing
with concurrent.futures.ProcessPoolExecutor() as executor:

## Multiprocessing to get alt kmer counts
    kmer_size = [3] * constrained_intervals.shape[0]
    ref_kmer_list = [ref_kmers] * constrained_intervals.shape[0]
    alt_kmer_list = [alt_kmer_dict] * constrained_intervals.shape[0]
    ccr_results = executor.map(process_constrained_intervals, constrained_intervals, ref_kmer_list, alt_kmer_list, kmer_size)
    
    ## Combine all the data
    my_list = []
    for result in ccr_results: 
        my_list.append(result)
        
    #my_list.head()
    ## Convert data into a pandas dataframe
    final = pd.DataFrame(my_list, columns = ['chromosome', 'start', 'stop', 'expected_mut_num', 'prob'])
    
    ## Write the data to an output file
    output_file = '/scratch/ucgd/lustre-work/quinlan/u1240855/somccr/data/output/constrained_intervals/constrained_processed_final.bed'
    final.to_csv(output_file, sep="\t", header=True, index=False)
    











