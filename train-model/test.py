#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 15:19:21 2021

@author: u1240855
"""
#%%
import pysam
import numpy as np
import json
import argparse 
import gzip 
import concurrent.futures
import functools
import copy
#%%

#%%
def get_regions_from_bed_file(regions_filename): 
  with gzip.open(regions_filename, mode='rt') as regions:
    return [bed_to_sam_string(region) for region in regions]
#%%

#%%
def compute_counts_concurrent():
  # https://docs.python.org/3/library/concurrent.futures.html#concurrent.futures.ProcessPoolExecutor
  with concurrent.futures.ProcessPoolExecutor() as executor: 
    regions = get_regions_from_bed_file(regions)
    return regions, list(executor.map(compute_counts_region, regions))
#%%

#%%
def compute_counts_region(region):  
  print_json({'region': region, **get_hostname_process_cpu()})

  args = parse_arguments()

  kmer_data = initialize_kmer_data(args) 
  
   # pysam.FastaFile uses the index produced by "samtools faidx": 
  # https://pysam.readthedocs.io/en/latest/api.html?highlight=fasta#pysam.FastaFile
  with pysam.TabixFile(args.mutations) as mutations, pysam.FastaFile(args.genome) as genome: 
    kmer_data = compute_total_counts_region(genome, kmer_data, region, args) 
    kmer_data = compute_snv_counts_region(mutations, genome, kmer_data, region, args)

  return kmer_data 
#%%

#%%
def compute_total_counts_region(genome, kmer_data, region, args): 
  # "fetch" API: https://pysam.readthedocs.io/en/latest/api.html?highlight=fasta#pysam.FastaFile
  # Note that fetch(region=region) does not work if the coordinates in "region" contains commas
  # Workaround is to parse "region" into "chromosome", "start", "end": 
  neutral_region = genome.fetch(*unpack(region))    

  print_string_as_info('Iterating over neutral region {} and counting kmers:'.format(region))
  number_of_kmers_containing_unspecified_bases = 0
  for position in np.arange(0, len(neutral_region), 1):
    try: 
      kmer = fetch_kmer_from_sequence(neutral_region, position, args.kmer_size)
    except IndexError:
      print_string_as_info_dim('IndexError at position: {}'.format(position))
      continue
    if contains_unspecified_bases(kmer): 
      number_of_kmers_containing_unspecified_bases += 1 
      continue 
    kmer_data[kmer]['cohort_sequence_count'] += args.number_samples 
    kmer_data[kmer]['sequence_count'] += 1

  print_string_as_info_dim(
    f'Number of kmers containing unspecified bases: {number_of_kmers_containing_unspecified_bases}'
  )
  print_unbuffered('')

  return kmer_data
#%%

#%%
def compute_snv_counts_region(mutations, genome, kmer_data, region, args):
  print_string_as_info('Fetching SNVs in region {} and incrementing corresponding kmer counts\n'.format(region))
  SNVs = fetch_SNVs(mutations, genome, region, args.__dict__)
  for SNV in SNVs: 
    kmer_data[SNV['kmer']]['ALT_counts'][SNV['ALT']] += 1
  return kmer_data
#%%

#%%
def estimate_mutation_probabilities_core(kmer_data):
  print_string_as_info('Estimating mutation probabilities\n')
  for kmer, data in kmer_data.items():
    probabilities = {}
    check_for_Ns(kmer) # sanity check 
    if data['sequence_count'] == 0: 
      for base in get_bases(): probabilities[base] = None
    else: 
      mutation_probability = 0.0
      for alternate_base in alternate_bases(kmer):
        # estimate probabilities for multinomial distribution: https://math.stackexchange.com/a/421838
        probabilities[alternate_base] = data['ALT_counts'][alternate_base]/data['cohort_sequence_count']
        mutation_probability += probabilities[alternate_base]
      probabilities[middle_base(kmer)] = 1.0 - mutation_probability
    data['estimated_mutation_probabilities'] = probabilities
    data['mutation_probability'] = mutation_probability
  return kmer_data
#%%

#%%
def combine_counts(x, y):
  result = copy.deepcopy(x)
  for kmer in result.keys(): 
    data_result, data_x, data_y = result[kmer], x[kmer], y[kmer]
    data_result['cohort_sequence_count'] = data_x['cohort_sequence_count'] + data_y['cohort_sequence_count']
    data_result['sequence_count'] = data_x['sequence_count'] + data_y['sequence_count']
    ALT_counts_result,  ALT_counts_x, ALT_counts_y = data_result['ALT_counts'], data_x['ALT_counts'], data_y['ALT_counts']
    for alternate_base in ALT_counts_result.keys(): 
        ALT_counts_result[alternate_base] = ALT_counts_x[alternate_base] + ALT_counts_y[alternate_base]
  return result
#%%

#%%
kmer_size = 3
genome = '/scratch/ucgd/lustre-work/quinlan/u1240855/constraint-tools/data/reference/grch38/hg38.analysisSet.fa.gz'
regions = '/scratch/ucgd/lustre-work/quinlan/u1240855/constraint-tools/dist/neutral-regions.bed.gz'
number_samples = 76156
mutations = '/scratch/ucgd/lustre-work/quinlan/u1240855/constraint-tools/data/gnomad/v3/gnomad_v3_variants.sorted.bed.gz'
model = '/scratch/ucgd/lustre-work/quinlan/u1240855/constraint-tools/tests'
training_mode = 'concurrent'
#%%

#%%

regions = get_regions_from_bed_file(regions)
genome = pysam.FastaFile(genome)
mutations = pysam.TabixFile(mutations)

## Get the kmer counts from the neutral region
kmer_data = compute_total_counts_region(genome, kmer_data, region, args) 

## Get kmer counts for variants wtihin each neutral region
SNVs = fetch_SNVs(mutations, genome, region, args.__dict__)


region = 'chr22:10510300-10510500'

kmer_data = functools.reduce(combine_counts, kmer_data)
kmer_data = estimate_mutation_probabilities_core(kmer_data) 
#%%