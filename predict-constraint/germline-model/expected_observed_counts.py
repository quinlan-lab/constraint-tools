import argparse 
import numpy as np 

# pysam API: 
# https://github.com/pysam-developers/pysam/blob/b82cbcae22c088e64fdb58f8acaf1e9773c7b088/pysam/libctabix.pyx
import pysam

from colorize import print_json  
from snvs import fetch_SNVs
from pack_unpack import unpack
from kmer import CpG, not_CpG
from windows import create_windows 
from get_p0s_p1s_p2s_p3s import get_p0s_p1s_p2s_p3s

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

# TODO: implement get_K_observed() 

# TODO: check 
def get_M(window, genome, mutations, model, number_chromosomes_min):
  return len(fetch_SNVs(mutations, genome, window['region'], meta=model, number_chromosomes_min=number_chromosomes_min))

def get_N_observed(window, genome, mutations, model):
  return len(fetch_SNVs(mutations, genome, window['region'], meta=model))

# https://github.com/quinlan-lab/constraint-tools/blob/main/define-model/germline-model.ipynb
def get_N_mean_variance_null(window, genome, model, log): 
  p0s_p1s_p2s_p3s = get_p0s_p1s_p2s_p3s(window, genome, model, log)

  mean_Ns = np.array([0*p0 + 1*p1 + 2*p2 + 3*p3 for p0, p1, p2, p3 in p0s_p1s_p2s_p3s])
  mean_N = np.sum(mean_Ns)

  mean_N2s = np.array([(0**2)*p0 + (1**2)*p1 + (2**2)*p2 + (3**2)*p3 for p0, p1, p2, p3 in p0s_p1s_p2s_p3s])
  variance_N = np.sum(mean_N2s - np.square(mean_Ns))

  return mean_N, variance_N

# TODO: 
def get_K_mean_variance_null(window, model, mutations, genome, number_chromosomes_min): 
  M_observed = get_M(window, genome, mutations, model, number_chromosomes_min)

  # RECALL TRAINING DOES THE FOLLOWING: 
  # SNVs = fetch_SNVs(mutations, genome, window['region'], args.__dict__, args.number_chromosomes_min)
  # SNV_count = len(SNVs)
  # singleton_count = len([SNV for SNV in SNVs if SNV['number_ALT_chromosomes'] == 1])
  # print_string_as_info_dim(f"{singleton_count} of {SNV_count} SNVs in {window['region']} are singletons")
  # singleton_counts[SNV_count][singleton_count] += 1

  print(model['singletonProbabilities'][M_observed]) 
  1/0
  return None, None

# TODO: 
def compute_Kbar_Kobserved(window, model, mutations, genome, number_chromosomes_min):
  K_mean_null, K_variance_null = get_K_mean_variance_null(window, model, mutations, genome, number_chromosomes_min)
  # IMPORTANT: report total number of counts (number of windows) used to estimate mean_K_null and variance_K_null 
  # Kbar = (K_observed - mean_K_null)/np.sqrt(variance_K_null)
  pass

def compute_Nbar_Nobserved(window, model, mutations, genome, log):
  N_mean_null, N_variance_null = get_N_mean_variance_null(window, genome, model, log)
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
    # TODO: finish: 
    # zip(*[compute_Kbar_Kobserved(window, model, mutations, genome) for window in windows])

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

if __name__ == '__main__':
  test()