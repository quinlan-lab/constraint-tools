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
import concurrent.futures
import collections

# using the "khmer" library might make this function a handful of times faster:
# https://khmer.readthedocs.io/en/v0.6.1-0/ktable.html
def compute_total_counts(genome, chromosome, start, end, kmer_data, args): 
  # "fetch" API: https://pysam.readthedocs.io/en/latest/api.html?highlight=fasta#pysam.FastaFile
  # Note that fetch(region=region) does not work if the coordinates in "region" contains commas
  # Workaround is to parse "region" into "chromosome", "start", "end": 
  
  print_string_as_info(f'Iterating over neutral region ({chromosome}:{start}-{end}) and counting reference kmers:')
  neutral_region = genome.fetch(chromosome, start, end)   
  neutral_sequence_length = len(neutral_region)
  for position in np.arange(0, neutral_sequence_length, 1):
    try: 
      kmer = fetch_kmer_from_sequence(neutral_region, position, args.kmer_size)
      if ('N' not in kmer) & ('M' not in kmer) & ('R' not in kmer): 
          #print(kmer)
          kmer_data[kmer]['cohort_sequence_count'] += args.number_tumors 
          kmer_data[kmer]['sequence_count'] += 1
    except IndexError:
      print_string_as_info_dim('IndexError at position: {}'.format(position))
      pass 
  print('')

  return kmer_data

def compute_snv_counts(mutations, genome, chromosome, start, end, kmer_data, args): 
  print_string_as_info('Fetching SNVs in region and incrementing corresponding kmer counts\n')
  SNVs = fetch_SNVs(mutations, genome, chromosome, start, end, args.__dict__)
  for SNV in SNVs: 
    kmer_data[SNV['kmer']]['ALT_counts'][SNV['ALT']] += 1
  return kmer_data   

def estimate_mutation_probabilities(neutral_region):
  args = parse_arguments()

  ## Initialize kmer data
  kmer_data = initialize_kmer_data(args) 
   
  '''
  Get the neutral region coordiantes
  '''
  ## Humanize variables
  chr = 0
  start = 1
  end = 2
  
  ## Get the chr, start, and end coordinates for the neutral region of interest
  chromosome = str(int(neutral_region[chr]))
  if chromosome == "23": 
      chromosome = "X"
  if chromosome == "24": 
      chromosome = "Y"
  start = int(neutral_region[start])
  end = int(neutral_region[end])
  
  # pysam.FastaFile uses the index produced by "samtools faidx": 
  # https://pysam.readthedocs.io/en/latest/api.html?highlight=fasta#pysam.FastaFile
  with pysam.TabixFile(args.mutations) as mutations, pysam.FastaFile(args.genome) as genome: 
    kmer_data = compute_total_counts(genome, chromosome, start, end, kmer_data, args)
    kmer_data = compute_snv_counts(mutations, genome, chromosome, start, end, kmer_data, args)  
  
  return kmer_data


def estimate_mutation_probabilities_core(kmer_data, args):
  print_string_as_info('Estimating mutation probabilities\n')
  for kmer, data in kmer_data.items():
    #print(f"printing kmer: {kmer}")
    #print(f"printing data: {data} ")
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

def parse_arguments(): 
  parser = argparse.ArgumentParser(description='')
  parser.add_argument('--kmer-size', type=int, dest='kmer_size', help='')
  parser.add_argument('--genome', type=str, help='')
  parser.add_argument('--neutral_regions', type=str, help='')
  parser.add_argument('--number-tumors', type=int, dest='number_tumors', help='')
  parser.add_argument('--output', type=str, help='')
  parser.add_argument('--mutations', type=str, help='')
  return parser.parse_args()


if __name__ == '__main__': 
  args = parse_arguments()
  
  ## Read in the neutral regions bed file 
  neutral_regions_file = args.neutral_regions
  neutral_regions = np.genfromtxt(neutral_regions_file, delimiter="\t")
  
  
  ## Multiprocessing - run kmer functions on each neutral region
  with concurrent.futures.ProcessPoolExecutor() as executor: 
      results = executor.map(estimate_mutation_probabilities, neutral_regions)
      
  '''
  Go through the results of multiprocessing to merge kmer counts from neutral regions...
  Combine alternate kmer counts calculated from the mutations file...

  Example output: 
  Neutral region 1, GAC kmer 
  {'AAA': {'CpG': False, 'cohort_sequence_count': 5850, 'sequence_count': 3, 'REF': 'A', 'ALT_counts': {'C': 0, 'G': 0, 'T': 0}}
   'ATA': {'CpG': False, 'cohort_sequence_count': 5850, 'sequence_count': 3, 'REF': 'A', 'ALT_counts': {'C': 0, 'G': 0, 'T': 0}}}
  
  Create nested dictionary: https://stackoverflow.com/questions/16333296/how-do-you-create-nested-dict-in-python 
  
  Increment python dictionary values with the same key/key combination: https://stackoverflow.com/questions/12992165/python-dictionary-increment
  '''
  
  print("Merging kmer data from each neutral region interval")
  
  kmer_dict = collections.defaultdict(dict)

  ## Iterate through kmer data from each neutral region interval 
  for result in results: 
      
      ## Iterate through each kmer 
      for key in result: 
          kmer = key
          kmer_dict[kmer]['CpG'] = result[kmer]['CpG']
          kmer_dict[kmer]['cohort_sequence_count'] = kmer_dict[kmer].get('cohort_sequence_count', 0) + result[kmer]['cohort_sequence_count']
          kmer_dict[kmer]['sequence_count'] = kmer_dict[kmer].get('sequence_count', 0) + result[kmer]['sequence_count']
          kmer_dict[kmer]['REF'] = result[kmer]['REF']   

          ## Iterate through each allele
          for key in result[kmer]['ALT_counts']: 
              allele = key
              if "ALT_counts" in kmer_dict[kmer]: 
                  kmer_dict[kmer]['ALT_counts'][allele] = kmer_dict[kmer]['ALT_counts'].get(allele, 0) + result[kmer]['ALT_counts'][allele]
              if "ALT_counts" not in kmer_dict[kmer]: 
                  kmer_dict[kmer]['ALT_counts'] = {str(allele): result[kmer]['ALT_counts'][allele]}
                  
         
  '''
  Estimate mutation probabilities based on reference and alternate kmer data
  '''  
  kmer_data = estimate_mutation_probabilities_core(kmer_dict, args)
  print(kmer_data)  
  
  '''
  Output data
  '''
  model_path = args.output + '/model.json'
  with open(model_path, 'w') as fh:
    json.dump({
      'mutations': args.mutations,
      'genome': args.genome,
      'kmer_size': args.kmer_size,
      'number_tumors': args.number_tumors,
      'kmer_data': kmer_data
    }, fh)
  print_string_as_info('Writing multinomial model to disk at: {}'.format(model_path))