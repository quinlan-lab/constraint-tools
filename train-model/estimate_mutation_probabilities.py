# pysam API: 
# https://github.com/pysam-developers/pysam/blob/b82cbcae22c088e64fdb58f8acaf1e9773c7b088/pysam/libctabix.pyx
import pysam

import numpy as np
import json
import argparse 

from kmer import compute_kmers, fetch_kmer_from_sequence
from colorize import print_json, print_string_as_info
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
      print_string_as_info('IndexError at position: {}'.format(position))
      pass 
  print('')

  return kmer_data   

def compute_snv_counts(mutations, genome, kmer_data, args): 
  print_string_as_info('Fetching SNVs in region and incrementing corresponding kmer counts\n')
  SNVs = fetch_SNVs(mutations, genome, args)
  for SNV in SNVs: 
    kmer_data[SNV['kmer']]['snv_count'] += 1
  return kmer_data   

def estimate_mutation_probabilities():
  args = parse_arguments()

  kmer_data = {
    kmer: {
      'cohort_sequence_count': 0,
      'sequence_count': 0,
      'snv_count': 0 
    } for kmer in compute_kmers(args.kmer_size)
  }

  # pysam.FastaFile uses the index produced by "samtools faidx": 
  # https://pysam.readthedocs.io/en/latest/api.html?highlight=fasta#pysam.FastaFile
  with pysam.TabixFile(args.mutations) as mutations, pysam.FastaFile(args.genome) as genome: 
    kmer_data = compute_total_counts(genome, kmer_data, args)
    kmer_data = compute_snv_counts(mutations, genome, kmer_data, args)    
    
  print_string_as_info('Estimating mutation probabilities\n')
  for kmer, data in kmer_data.items():
      try: 
          kmer_data[kmer]['estimated_mutation_probability'] = data['snv_count']/data['cohort_sequence_count']
      except ZeroDivisionError: 
          kmer_data[kmer]['estimated_mutation_probability'] = None 
          print_string_as_info('ZeroDivisionError:')
          print_string_as_info(kmer)
          print_json(data)
      kmer_data[kmer]['number_tumors'] = args.number_tumors
      
  model_path = args.output + '/binomial_model.json'
  with open(model_path, 'w') as fh:
    json.dump(kmer_data, fh)
  print_string_as_info('Writing binomial model to disk at: {}'.format(model_path))

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