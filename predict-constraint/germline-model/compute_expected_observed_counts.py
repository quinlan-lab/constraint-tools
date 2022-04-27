import argparse 
import collections
import numpy as np 

# pysam API: 
# https://github.com/pysam-developers/pysam/blob/b82cbcae22c088e64fdb58f8acaf1e9773c7b088/pysam/libctabix.pyx
import pysam

from colorize import print_json  
from snvs import fetch_SNVs
from pack_unpack import unpack
from kmer import (
  CpG, 
  not_CpG, 
  fetch_kmers,
  compute_possible_ALT_states_core
)
from windows import create_windows 

import color_traceback

def compute_SNV_positions_frequencies_filtered(SNVs, filter_function):
  SNV_positions_frequencies = [
    (SNV['position'], SNV['number_ALT_chromosomes']) 
    for SNV in SNVs if filter_function(SNV['kmer'])
  ]  
  return tuple(zip(*SNV_positions_frequencies))

def compute_SNV_positions_frequencies(mutations, genome, region, model, number_chromosomes_min):
  SNVs = fetch_SNVs(mutations, genome, region, meta=model, number_chromosomes_min=number_chromosomes_min)
  return (
    compute_SNV_positions_frequencies_filtered(SNVs, CpG),
    compute_SNV_positions_frequencies_filtered(SNVs, not_CpG)
  )

def pull_element(list_, index): 
  try:
    return list_[index]
  except IndexError: 
    return []

# TODO: 
def fetch_distribution_K(window, genome, model):
  # 1. compute observed SNV count for given window 
  # 2. pull out the singleton-count distribution for that SNV count (xs, ys)
  # 3. return xs, ys, and counts (where counts is the number of windows used to compute the distribution)
  pass

def create_column_of_Ns_core(number_examples, p0_p1_p2_p3): 
  return np.random.choice(a=[0, 1, 2, 3], size=(number_examples, 1), p=p0_p1_p2_p3)

def fetch_distribution_N(window, genome, model): 
  number_examples = 100000
  create_column_of_Ns = lambda p0_p1_p2_p3: create_column_of_Ns_core(number_examples, p0_p1_p2_p3)
  list_of_columns_of_Ns = map(create_column_of_Ns, get_p0s_p1s_p2s_p3s(window, genome, model))
  Ns = np.column_stack(list(list_of_columns_of_Ns))
  N = np.sum(Ns, axis=1)
  N_histogram = collections.Counter(N)    
  max_N = max(N)
  xs = range(max_N+1)
  ys = np.array([
      N_histogram[x] 
      if x in N_histogram.keys()
      else 0 
      for x in xs
  ])/number_examples  
  return {
    'n': list(xs), # list makes range JSON-serializable
    'p(n)': list(ys) # list makes numpy array JSON-serializable
  }

# TODO: 
def compute_Kbar_Kobserved(): 
  # compute mean_K_null 
  # compute variance_K_null 
  # IMPORTANT: report total number of counts (number of windows) used to estimate mean_K_null and variance_K_null 
  # Kbar = (K_observed - mean_K_null)/np.sqrt(variance_K_null)
  pass

def get_N_observed(window, genome, mutations, model):
  return len(fetch_SNVs(mutations, genome, window['region'], meta=model))

def get_substitution_probability(model, kmer, ALT_multiplicity): 
  return sum(
    model['kmerProbabilities'][kmer][alt_state] 
    for alt_state in compute_possible_ALT_states_core(kmer, ALT_multiplicity)
  )

# https://github.com/quinlan-lab/constraint-tools/blob/main/define-model/germline-model.ipynb
def get_p0s_p1s_p2s_p3s(window, genome, model, log=True):
  p0s_p1s_p2s_p3s = []
  for kmer in fetch_kmers(window['region'], genome, model['kmerSize'], log): 
    p1_p2_p3 = [get_substitution_probability(model, kmer, ALT_multiplicity) for ALT_multiplicity in [1, 2, 3]]
    p0 = 1 - sum(p1_p2_p3) 
    p0s_p1s_p2s_p3s.append([p0] + p1_p2_p3)
  return p0s_p1s_p2s_p3s

# https://github.com/quinlan-lab/constraint-tools/blob/main/define-model/germline-model.ipynb
def get_null_mean_variance(window, genome, model, log): 
  p0s_p1s_p2s_p3s = get_p0s_p1s_p2s_p3s(window, genome, model, log)

  mean_Ns = np.array([0*p0 + 1*p1 + 2*p2 + 3*p3 for p0, p1, p2, p3 in p0s_p1s_p2s_p3s])
  mean_N = np.sum(mean_Ns)

  mean_N2s = np.array([(0**2)*p0 + (1**2)*p1 + (2**2)*p2 + (3**2)*p3 for p0, p1, p2, p3 in p0s_p1s_p2s_p3s])
  variance_N = np.sum(mean_N2s - np.square(mean_Ns))

  return mean_N, variance_N

def compute_Nbar_Nobserved(window, model, mutations, genome, log):
  N_mean_null, N_variance_null = get_null_mean_variance(window, genome, model, log)
  N_observed = get_N_observed(window, genome, mutations, model)
  N_bar = (N_observed - N_mean_null)/np.sqrt(N_variance_null)
  return N_bar, N_observed

def compute_expected_observed_counts(region, model, window_stride, number_chromosomes_min, log=True):
  with pysam.TabixFile(model['mutations']) as mutations, pysam.FastaFile(model['genome']) as genome:
    windows = create_windows(model['windowSize'], window_stride, region, genome, region_contains_windows=True)    
    N_bars, N_observeds = zip(*[compute_Nbar_Nobserved(window, model, mutations, genome, log) for window in windows])
    (
      SNV_positions_frequencies_CpG_positive, 
      SNV_positions_frequencies_CpG_negative
     ) = compute_SNV_positions_frequencies(mutations, genome, region, model, number_chromosomes_min)

  chromosome, start, end = unpack(region)

  return {
    'region': region,
    'chromosome': chromosome,
    'start': start,
    'end': end,
    'windowPositions': [window['position'] for window in windows], 
    'windowRegions': [window['region'] for window in windows], 
    'NBars': N_bars, 
    'NObserveds': N_observeds, 
    'snvCpGPositivePositions': pull_element(SNV_positions_frequencies_CpG_positive, index=0), 
    'snvCpGPositiveFrequencies': pull_element(SNV_positions_frequencies_CpG_positive, index=1), 
    'snvCpGNegativePositions': pull_element(SNV_positions_frequencies_CpG_negative, index=0),
    'snvCpGNegativeFrequencies': pull_element(SNV_positions_frequencies_CpG_negative, index=1)
  }

def parse_arguments(): 
  parser = argparse.ArgumentParser(description='')
  parser.add_argument('--region', type=str, help='')
  parser.add_argument('--model', type=str, help='')
  parser.add_argument('--window-stride', type=int, help='', dest='window_stride')
  return parser.parse_args()

def test(): 
  args = parse_arguments()   

  from read_model import read_model
  model = read_model(args.model)

  print_json(compute_expected_observed_counts(
    args.region, 
    model, 
    args.window_stride, 
    number_chromosomes_min=model['numberChromosomesMin']
  ))

  with pysam.FastaFile(model['genome']) as genome:
    windows = create_windows(model['windowSize'], args.window_stride, args.region, genome, region_contains_windows=True)      
    window = windows[0]
    print_json(fetch_distribution_N(window, genome, model))

if __name__ == '__main__':
  test()