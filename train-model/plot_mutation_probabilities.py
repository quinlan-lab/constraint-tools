import json
import argparse 
import numpy as np 

from kmer import initialize_kmer_data, fetch_kmer_from_sequence, alternate_bases, middle_base
from colorize import print_json, print_string_as_info, print_string_as_info_dim
import color_traceback 
from fetch_SNVs import fetch_SNVs 

def parse(region): 
  chromosome, start_end = region.split(':')
  start, end = map(lambda s: int(s.replace(',', '')), start_end.split('-'))
  return chromosome, start, end

def compute_total_counts(genome, kmer_data, args): 
  # "fetch" API: https://pysam.readthedocs.io/en/latest/api.html?highlight=fasta#pysam.FastaFile
  # Note that fetch(region=region) does not work if the coordinates in "region" contains commas
  # Workaround is to parse "region" into "chromosome", "start", "end": 
  neutral_region = genome.fetch(*parse(args.region)).upper()    

  print_string_as_info('Iterating over neutral region and counting kmers:')
  for position in np.arange(0, len(neutral_region), 1):
    try: 
      kmer = fetch_kmer_from_sequence(neutral_region, position, args.kmer_size)
      kmer_data[kmer]['cohort_sequence_count'] += args.number_tumors 
      kmer_data[kmer]['sequence_count'] += 1
    except IndexError:
      print_string_as_info_dim('IndexError at position: {}'.format(position))
      pass 
  print('')

  return kmer_data   

def compute_snv_counts(mutations, genome, kmer_data, args): 
  print_string_as_info('Fetching SNVs in region and incrementing corresponding kmer counts\n')
  SNVs = fetch_SNVs(mutations, genome, args)
  for SNV in SNVs: 
    kmer_data[SNV['kmer']]['ALT_counts'][SNV['ALT']] += 1
  return kmer_data   

def estimate_mutation_probabilities_core(kmer_data, args):
  print_string_as_info('Estimating mutation probabilities\n')
  for kmer, data in kmer_data.items():
    probabilities = {}
    for alternate_base in alternate_bases(kmer):
      try: 
        probabilities[alternate_base] = data['ALT_counts'][alternate_base]/data['cohort_sequence_count']
      except ZeroDivisionError: 
        probabilities[alternate_base] = None 
        print_string_as_info('ZeroDivisionError:')
        print_string_as_info(kmer)
        print_json(data)
    probabilities[middle_base(kmer)] = 1.0 - np.sum(
      [probabilities[alternate_base] for alternate_base in alternate_bases(kmer)]
    )
    data['estimated_mutation_probabilities'] = probabilities
    data['number_tumors'] = args.number_tumors
  return kmer_data

def estimate_mutation_probabilities():
  args = parse_arguments()

  kmer_data = initialize_kmer_data(args) 

  # pysam.FastaFile uses the index produced by "samtools faidx": 
  # https://pysam.readthedocs.io/en/latest/api.html?highlight=fasta#pysam.FastaFile
  with pysam.TabixFile(args.mutations) as mutations, pysam.FastaFile(args.genome) as genome: 
    kmer_data = compute_total_counts(genome, kmer_data, args)
    kmer_data = compute_snv_counts(mutations, genome, kmer_data, args)    
    
  kmer_data = estimate_mutation_probabilities_core(kmer_data, args) 

  model_path = args.output + '/multinomial_model.json'
  with open(model_path, 'w') as fh:
    json.dump(kmer_data, fh)
  print_string_as_info('Writing multinomial model to disk at: {}'.format(model_path))

