# pysam API: 
# https://github.com/pysam-developers/pysam/blob/b82cbcae22c088e64fdb58f8acaf1e9773c7b088/pysam/libctabix.pyx
import pysam
import numpy as np
import argparse 
from kmer import CpG, not_CpG, fetch_kmer_from_genome, compute_left_right
import json
from colorize import print_json, print_string_as_info, print_string_as_info_dim


'''
Version 1: 
    Function to calculate probability of 0 mutations within specific region of putative purifying selection
    Step 1: calculate expected mutation counts within the RoI based on the RoI's kmer distribution
    Step 2: calculate observed mutation counts within the RoI based on MAF mutation file
    Step 3: use the expected mutation count to get probability of purifying selection (n=0 mutations) within the RoI
'''

## Compute expected mutation counts within the region of interest --> use binomial statistics
def compute_window_expected_mutation_count(genome, args): 
    region = genome.fetch(chromosome, start, end)
    region_seq_length = len(region)
    print(region)
    print(region_seq_length)

    # for position in np.arrange(0, region_seq_length, 1): 
    #     try: 
    #         kmer = fetch_kmer_from_sequene(args.chromosome, args.start, args.end, position, args.kmer_size)
    #     except IndexError:
    #         print_string_as_info_dim('IndexError at position:{}'.format(position))
    #         pass
    # print('')

    # print(kmer)
    # return kmer

## Compute observed mutation counts within the region of interest
#def compute_window_observed_mutation_count(region, genome, kmer_size):
    

## Calculate probability of purifying selection within the region of interest
#def compute_probability(expected_count, observed_count):
    

## Define input arguments
def parse_arguments(): 
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--kmer-size', type=int, dest='kmer_size', help='')
    parser.add_argument('--genome', type=str, help='')
    parser.add_argument('--chromosome', type=str, help='')
    parser.add_argument('--start', type=str, help='')
    parser.add_argument('--end', type=str, help='')
    parser.add_argument('--model', type=str, help='')
    parser.add_argument('--number-tumors', type=int, dest='number_tumors', help='')
    parser.add_argument('--output', type=str, help='')
    parser.add_argument('--mutations', type=str, help='')
    return parser.parse_args()

## Define main function
if __name__ == '__main__': 
    args = parse_arguments()
    
    with pysam.TabixFile(args.mutations) as mutations, pysam.FastaFile(args.genome) as genome, json.load(model) as model: 
        compute_window_expected_mutation_count(genome, model, args)
      

  

  




  

  
