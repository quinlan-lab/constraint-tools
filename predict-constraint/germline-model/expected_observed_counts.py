import argparse 
import numpy as np 

# pysam API: 
# https://github.com/pysam-developers/pysam/blob/b82cbcae22c088e64fdb58f8acaf1e9773c7b088/pysam/libctabix.pyx
import pysam

import pyranges as pr

from colorize import print_json  
from snvs import fetch_SNVs
from pack_unpack import unpack
from kmer import CpG, not_CpG
from windows import create_windows 
from get_p0s_p1s_p2s_p3s import get_p0s_p1s_p2s_p3s
from get_M_K import get_M_K_testTime 
from neutral_regions import get_all_neutral_regions
from exons import get_canonical_exons

import color_traceback

def compute_SNV_positions_frequencies_filtered(SNVs, filter_function):
  SNV_positions_frequencies = [
    (SNV['position'], SNV['number_ALT_chromosomes']) 
    for SNV in SNVs if filter_function(SNV['kmer'])
  ]  
  return tuple(zip(*SNV_positions_frequencies))

def compute_SNV_positions_frequencies(mutations, genome, region, model):
  SNVs = fetch_SNVs(mutations, genome, region, meta=model, number_chromosomes_min=model['numberChromosomesMin'])
  return (
    compute_SNV_positions_frequencies_filtered(SNVs, CpG),
    compute_SNV_positions_frequencies_filtered(SNVs, not_CpG)
  )

def pull_element(list_, index): 
  try:
    return list_[index]
  except IndexError: 
    return []

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

def get_K_mean_variance_null(M, model): 
  possible_Ks = np.arange(0, M+1)
  window_counts = model['singletonCounts'][M]

  # https://stackoverflow.com/a/50786849/6674256
  K_mean_null = np.average(possible_Ks, weights=window_counts, returned=False)
  K_variance_null = np.average((possible_Ks - K_mean_null)**2, weights=window_counts, returned=False)

  return K_mean_null, K_variance_null

def compute_Kbar_Kobserved_M(window, model, mutations, genome):
  M, K_observed = get_M_K_testTime(window, mutations, genome, model)
  if M > 0 and M in model['singletonCounts'].keys(): 
    K_mean_null, K_variance_null = get_K_mean_variance_null(M, model)
    K_bar = (K_observed - K_mean_null)/np.sqrt(K_variance_null)
  else: 
    K_bar = None 
  return K_bar, K_observed, M

def compute_Nbar_Nobserved(window, model, mutations, genome, log):
  N_mean_null, N_variance_null = get_N_mean_variance_null(window, genome, model, log)
  N_observed = get_N_observed(window, genome, mutations, model)
  N_bar = (N_observed - N_mean_null)/np.sqrt(N_variance_null)
  return N_bar, N_observed

def filter_by_regions(windows, N_bars, K_bars, regions, how): 
  chromosomes, starts, ends = tuple(zip(*[unpack(window['region']) for window in windows]))

  # https://biocore-ntnu.github.io/pyranges/loadingcreating-pyranges.html
  z_scores = pr.PyRanges(chromosomes=chromosomes, starts=starts, ends=ends)

  # https://biocore-ntnu.github.io/pyranges/manipulating-the-data-in-pyranges.html
  z_scores.windowPositions = [window['position'] for window in windows]
  z_scores.NBars = N_bars
  z_scores.KBars = K_bars
  
  # https://biocore-ntnu.github.io/pyranges/intersecting-ranges.html
  z_scores_filtered = z_scores.intersect(regions, how=how)

  if z_scores_filtered.df.empty: 
    return None, None, None

  z_scores_filtered.KBars = z_scores_filtered.KBars.replace({np.nan: None})

  return (
    z_scores_filtered.windowPositions.tolist(),
    z_scores_filtered.NBars.tolist(), 
    z_scores_filtered.KBars.tolist()
  )

def compute_expected_observed_counts(region, model, window_stride, log=True):
  with pysam.TabixFile(model['mutations']) as mutations, pysam.FastaFile(model['genome']) as genome:
    windows = create_windows(model['windowSize'], window_stride, region, genome, region_contains_windows=True)
    N_bars, N_observeds = zip(*[compute_Nbar_Nobserved(window, model, mutations, genome, log) for window in windows])
    (
      SNV_positions_frequencies_CpG_positive, 
      SNV_positions_frequencies_CpG_negative
     ) = compute_SNV_positions_frequencies(mutations, genome, region, model)
    K_bars, K_observeds, Ms = zip(*[
      compute_Kbar_Kobserved_M(window, model, mutations, genome) 
      for window in windows
    ])

  chromosome, start, end = unpack(region)

  (
    window_positions_neutral_regions, 
    N_bars_neutral_regions, 
    K_bars_neutral_regions
  ) = filter_by_regions(windows, N_bars, K_bars, regions=get_all_neutral_regions(model), how='containment')
  (
    window_positions_exons, 
    N_bars_exons, 
    K_bars_exons
  ) = filter_by_regions(windows, N_bars, K_bars, regions=get_canonical_exons(), how=None)

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
    'snvCpGNegativeFrequencies': pull_element(SNV_positions_frequencies_CpG_negative, index=1),
    'KBars': K_bars, 
    'KObserveds': K_observeds, 
    'Ms': Ms,
    'windowPositionsNeutralRegions': window_positions_neutral_regions,
    'NBarsNeutralRegions': N_bars_neutral_regions,
    'KBarsNeutralRegions': K_bars_neutral_regions,
    'windowPositionsExons': window_positions_exons,
    'NBarsExons': N_bars_exons,
    'KBarsExons': K_bars_exons,
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
    args.window_stride
  ))

if __name__ == '__main__':
  test()