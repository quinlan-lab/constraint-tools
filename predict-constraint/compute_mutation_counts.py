# pysam API: 
# https://github.com/pysam-developers/pysam/blob/b82cbcae22c088e64fdb58f8acaf1e9773c7b088/pysam/libctabix.pyx
import pysam

import json 
from scipy.stats import binom
import collections
import numpy as np 

from colorize import print_string_as_info_dim
from snvs import fetch_SNVs
from pack_unpack import unpack
from kmer import CpG, not_CpG, fetch_kmer_from_genome
from windows import create_windows 

import color_traceback

def compute_expected_mutation_count(kmer, model): 
  return binom.stats(
    model['number_samples'], 
    model['kmer_data'][kmer]['mutation_probability'], 
    moments = 'm'
  )

def compute_window_expected_mutation_count(window, model, genome): 
  chromosome, start, end = unpack(window['region'])
  window_expected_mutation_count = 0
  for position in np.arange(start, end, 1):
    try: 
      kmer = fetch_kmer_from_genome(genome, chromosome, position, model['kmer_size'])
      window_expected_mutation_count += compute_expected_mutation_count(kmer, model)
    except ValueError:
      print_string_as_info_dim('ValueError at position {} in window {}'.format(position, window))
  return window_expected_mutation_count

# # import tensorflow.compat.v2 as tf
# import tensorflow_probability as tfp
# # tf.enable_v2_behavior()
# tfd = tfp.distributions 
# from collections import Counter

# def estimate_probability_of_number_of_mutations_in_interval(
#   number_kmers,
#   kmer_mutation_probabilities,
#   number_examples = 1000 # number of samples to draw from each binomial distribution 
# ): 
#   number_of_mutations_of_kmers = tfd.Binomial(
#     total_count = number_kmers, 
#     probs = kmer_mutation_probabilities
#   ).sample(number_examples)
  
#   number_of_mutations_in_interval = tf.reduce_sum(number_of_mutations_of_kmers, axis=1).numpy()

#   counts_of_number_of_mutations_in_interval = Counter(number_of_mutations_in_interval)    
#   max_number_of_mutations_in_interval = int(np.max(number_of_mutations_in_interval))
#   probability_of_number_of_mutations_in_interval = np.array([
#     counts_of_number_of_mutations_in_interval[value] 
#     if value in counts_of_number_of_mutations_in_interval 
#     else 0 
#     for value in range(max_number_of_mutations_in_interval)
#   ])/number_examples
  
#   return max_number_of_mutations_in_interval, probability_of_number_of_mutations_in_interval

# def compute_window_expected_mutation_count_distribution(window, model, genome): 
#   chromosome, start, end = unpack(window['region'])
#   number_kmers = []
#   kmer_mutation_probabilities = []
#   for position in np.arange(start, end, 1):
#     try: 
#       kmer = fetch_kmer_from_genome(genome, chromosome, position, model['kmer_size'])
#       number_kmers.append(model['number_samples'])
#       kmer_mutation_probabilities.append(model['kmer_data'][kmer]['mutation_probability'])
#     except ValueError:
#       print_string_as_info_dim('ValueError at position {} in window {}'.format(position, window))

#   max_number_of_mutations_in_interval, probability_of_number_of_mutations_in_interval = \
#     estimate_probability_of_number_of_mutations_in_interval(number_kmers, kmer_mutation_probabilities)
  
#   print_string_as_info_dim(max_number_of_mutations_in_interval)
#   print_string_as_info_dim(probability_of_number_of_mutations_in_interval)
#   1/0

#   return max_number_of_mutations_in_interval, probability_of_number_of_mutations_in_interval

def compute_window_expected_mutation_counts(windows, model, genome): 
  return [compute_window_expected_mutation_count(window, model, genome) for window in windows]
  
# def compute_window_expected_mutation_count_distributions(windows, model, genome): 
#   return [compute_window_expected_mutation_count_distribution(window, model, genome) for window in windows]
  
def compute_observed_mutation_count(mutations, genome, region, model):
  return len(fetch_SNVs(mutations, genome, region, meta=model))

def compute_window_observed_mutation_counts(mutations, genome, windows, model):
  return [compute_observed_mutation_count(mutations, genome, window['region'], model) for window in windows]

def compute_SNV_positions_frequencies(SNVs, filter_function):
  SNV_positions = [SNV['position'] for SNV in SNVs if filter_function(SNV['kmer'])]  
  # find frequency of each unique coordinate using Counter: 
  # https://stackoverflow.com/a/2162045/6674256  
  SNV_positions_frequencies = collections.Counter(SNV_positions)
  return tuple(zip(*SNV_positions_frequencies.items()))

def compute_lollipops(mutations, genome, region, model):
  SNVs = fetch_SNVs(mutations, genome, region, meta=model)
  return (
    compute_SNV_positions_frequencies(SNVs, CpG),
    compute_SNV_positions_frequencies(SNVs, not_CpG)
  )

def pull_element(list_, index): 
  try:
    return list_[index]
  except IndexError: 
    return []

# TODO: make this another flask endpoint that is hit when one wants to plot the distribution of singleton count in the web app: 
  # plt.plot(x_K, probability_of_K, marker='s')
  # plt.title(f'number of ALT alleles in interval = {m}')
  # plt.xlabel(f'number of singletons in interval, $k$')
  # _ = plt.ylabel('probability, $P[K=k]$')
def fetch_distribution_K(m): # m = number of ALT alleles in interval 
  return x_K, probability_of_K

def compute_p0s_p1s_p2s_p3s(): 
  # TODO: 
  # for each site, do the following: 
  # if site has a kmer with an unspecified base, then ignore site
  # otherwise, pull out the kmer polymorphism probabilities, and use them to compute p0, p1, etc
  # e.g., p1 = np.sum([kmer_probabilities[kmer][ALT_state] for ALT_state in compute_ALT_states(kmer, ALT_multiplicity=1)])
  # p0, p1, etc are defined in the section entitled "A model to predict the number of ALT alleles in a genomic interval" at https://github.com/quinlan-lab/constraint-tools/blob/main/define-model/germline-model.ipynb

# TODO: make this another flask endpoint that is hit when one wants to plot the distribution of SNV count in the web app: 
#   plt.plot(x_N, probability_of_N, marker='s')
#  plt.title(f'interval length = {number_sites}')
#  plt.xlabel(f'number of ALT alleles in interval, $n$')
#  _ = plt.ylabel('probability, $P[N=n]$')
def fetch_distribution_N(): 
  import numpy as np
  from collections import Counter

  N = []
  number_examples = 10000
  for example in range(number_examples):
    Ns = [np.random.choice(a=[0, 1, 2, 3], p=p0_p1_p2_p3) for p0_p1_p2_p3 in compute_p0s_p1s_p2s_p3s()]
    N.append(np.sum(Ns))

  N_histogram = Counter(N)    
  max_N = max(N)
  x_N = range(max_N+1)
  probability_of_N = np.array([
      N_histogram[value] 
      if value in N_histogram 
      else 0 
      for value in x_N
  ])/number_examples
  
  return x_N, probability_of_N

def compute_Nbar(): 
  # TODO: 
  # compute N_observed = number of polymorphic sites observed in a window
  # compute mean of N under null model: 
  # mean_Ns = np.array([0*p0 + 1*p1 + 2*p2 + 3*p3 for p0, p1, p2, p3 in p0s_p1s_p2s_p3s])
  # mean_N_null = np.sum(mean_Ns)
  # compute variance under null:
  # mean_N2s = np.array([(0**2)*p0 + (1**2)*p1 + (2**2)*p2 + (3**2)*p3 for p0, p1, p2, p3 in p0s_p1s_p2s_p3s])
  # variance_N_null = np.sum(mean_N2s - np.square(mean_Ns))
  # compute Nbar: 
  # N_observed = np.array(N_observed)
  # Nbar = (N_observed - mean_N_null)/np.sqrt(variance_N_null)

# TODO: 
def compute_Kbar(): 
  # compute mean_K_null 
  # compute variance_K_null 
  # Kbar = (K_observed - mean_K_null)/np.sqrt(variance_K_null)


def compute_mutation_counts(region, model_filename, window_size, window_stride):
  with open(model_filename) as fh:
    model = json.load(fh)

  with pysam.TabixFile(model['mutations']) as mutations, pysam.FastaFile(model['genome']) as genome:
    windows = create_windows(window_size, window_stride, region, genome)
    window_expected_mutation_counts = compute_window_expected_mutation_counts(windows, model, genome)
    # window_expected_mutation_count_distributions = compute_window_expected_mutation_count_distributions(windows, model, genome)
    window_observed_mutation_counts = compute_window_observed_mutation_counts(mutations, genome, windows, model)
    lollipops_CpG_positive, lollipops_CpG_negative = compute_lollipops(mutations, genome, region, model)

  chromosome, start, end = unpack(region)

  return {
    'region': region,
    'chromosome': chromosome,
    'start': start,
    'end': end,
    'windowPositions': [window['position'] for window in windows], 
    'windowRegions': [window['region'] for window in windows], 
    'windowExpectedMutationCounts': window_expected_mutation_counts, 
    'windowObservedMutationCounts': window_observed_mutation_counts, 
    'lollipopsCpGPositivePositions': pull_element(lollipops_CpG_positive, index=0), 
    'lollipopsCpGPositiveHeights': pull_element(lollipops_CpG_positive, index=1), 
    'lollipopsCpGNegativePositions': pull_element(lollipops_CpG_negative, index=0),
    'lollipopsCpGNegativeHeights': pull_element(lollipops_CpG_negative, index=1)
  }

import argparse 

def parse_arguments(): 
  parser = argparse.ArgumentParser(description='')
  parser.add_argument('--region', type=str, help='')
  parser.add_argument('--model', type=str, help='')
  parser.add_argument('--window-size', type=int, help='', dest='window_size')
  parser.add_argument('--window-stride', type=int, help='', dest='window_stride')
  return parser.parse_args()

if __name__ == '__main__':
  args = parse_arguments() 
  
  print(
    json.dumps(
      compute_mutation_counts(
        region=args.region, 
        model_filename=args.model, 
        window_size=args.window_size, 
        window_stride=args.window_stride
      ), 
      indent=2
    )
  )