# pysam API: 
# https://github.com/pysam-developers/pysam/blob/b82cbcae22c088e64fdb58f8acaf1e9773c7b088/pysam/libctabix.pyx
import pysam

import numpy as np
import json
import argparse 

from kmer import initialize_kmer_data, fetch_kmer_from_sequence, alternate_bases, middle_base, get_bases
from colorize import print_json, print_string_as_info, print_string_as_info_dim
import color_traceback 
from fetch_SNVs import fetch_SNVs 
from pack_unpack import unpack 

# using the "khmer" library might make this function a handful of times faster:
# https://khmer.readthedocs.io/en/v0.6.1-0/ktable.html
def compute_total_counts(genome, kmer_data, args): 
  # "fetch" API: https://pysam.readthedocs.io/en/latest/api.html?highlight=fasta#pysam.FastaFile
  # Note that fetch(region=region) does not work if the coordinates in "region" contains commas
  # Workaround is to parse "region" into "chromosome", "start", "end": 
  neutral_region = genome.fetch(*unpack(args.region))    
  neutral_sequence_length = len(neutral_region)

  print_string_as_info('Iterating over neutral region and counting kmers:')
  for position in np.arange(0, neutral_sequence_length, 1):
    try: 
      kmer = fetch_kmer_from_sequence(neutral_region, position, args.kmer_size)
      kmer_data[kmer]['cohort_sequence_count'] += args.number_tumors 
      kmer_data[kmer]['sequence_count'] += 1
    except IndexError:
      print_string_as_info_dim('IndexError at position: {}'.format(position))
      pass 
  print('')

  return kmer_data, neutral_sequence_length

def compute_snv_counts(mutations, genome, kmer_data, args): 
  print_string_as_info('Fetching SNVs in region and incrementing corresponding kmer counts\n')
  SNVs = fetch_SNVs(mutations, genome, args.region, args.__dict__)
  for SNV in SNVs: 
    kmer_data[SNV['kmer']]['ALT_counts'][SNV['ALT']] += 1
  return kmer_data   

def estimate_mutation_probabilities_core(kmer_data, args):
  print_string_as_info('Estimating mutation probabilities\n')
  for kmer, data in kmer_data.items():
    probabilities = {}
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

def estimate_mutation_probabilities():
  args = parse_arguments()

  kmer_data = initialize_kmer_data(args) 

  # pysam.FastaFile uses the index produced by "samtools faidx": 
  # https://pysam.readthedocs.io/en/latest/api.html?highlight=fasta#pysam.FastaFile
  with pysam.TabixFile(args.mutations) as mutations, pysam.FastaFile(args.genome) as genome: 
    kmer_data, neutral_sequence_length = compute_total_counts(genome, kmer_data, args)
    kmer_data = compute_snv_counts(mutations, genome, kmer_data, args)    
    
  kmer_data = estimate_mutation_probabilities_core(kmer_data, args) 

  model_path = args.output + '/model.json'
  with open(model_path, 'w') as fh:
    json.dump({
      'mutations': args.mutations,
      'genome': args.genome,
      'kmer_size': args.kmer_size,
      'number_tumors': args.number_tumors,
      'neutral_sequence_length': neutral_sequence_length,
      'kmer_data': kmer_data
    }, fh)
  print_string_as_info('Writing multinomial model to disk at: {}'.format(model_path))

def parse_arguments(): 
  parser = argparse.ArgumentParser(description='')
  parser.add_argument('--kmer-size', type=int, dest='kmer_size', help='')
  parser.add_argument('--genome', type=str, help='')
  parser.add_argument('--region', type=str, help='')
  parser.add_argument('--number-tumors', type=int, dest='number_tumors', help='')
  parser.add_argument('--output', type=str, help='')
  parser.add_argument('--mutations', type=str, help='')
  return parser.parse_args()

if __name__ == '__main__': 
  estimate_mutation_probabilities()  