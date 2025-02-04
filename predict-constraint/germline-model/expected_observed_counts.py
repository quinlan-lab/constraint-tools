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
from get_M_K import get_M_K_testTime 
from trustworthy_noncoding_regions import get_all_trustworthy_noncoding_regions
from exons import get_canonical_exons
from compute_Nbar_Nobserved import compute_Nbar_Nobserved 

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

def get_K_mean_variance_null(M, model): 
  possible_Ks = np.arange(0, M+1)
  window_counts = model['singletonCounts'][M]

  # https://stackoverflow.com/a/50786849/6674256
  K_mean_null = np.average(possible_Ks, weights=window_counts, returned=False)
  K_variance_null = np.average((possible_Ks - K_mean_null)**2, weights=window_counts, returned=False)

  return K_mean_null, K_variance_null

def compute_Kbar_Kobserved_M(window, model, mutations, genome):
  M, K_observed = get_M_K_testTime(window, mutations, genome, model)

  if M == 0: 
    return None, K_observed, M
  if M not in model['singletonCounts'].keys(): 
    return None, K_observed, M
  number_null_windows = np.sum(model['singletonCounts'][M])
  if number_null_windows < 50: 
    return None, K_observed, M

  K_mean_null, K_variance_null = get_K_mean_variance_null(M, model)
  K_bar = (K_observed - K_mean_null)/np.sqrt(K_variance_null)
  return K_bar, K_observed, M

def filter_by_regions(windows, N_bars, N_observeds, K_bars, K_observeds, regions, how): 
  chromosomes, starts, ends = tuple(zip(*[unpack(window['region']) for window in windows]))

  # https://biocore-ntnu.github.io/pyranges/loadingcreating-pyranges.html
  z_scores = pr.PyRanges(chromosomes=chromosomes, starts=starts, ends=ends)

  # https://biocore-ntnu.github.io/pyranges/manipulating-the-data-in-pyranges.html
  z_scores.windowPositions = [window['position'] for window in windows]
  z_scores.NBars = N_bars
  z_scores.NObserveds = N_observeds
  z_scores.KBars = K_bars
  z_scores.KObserveds = K_observeds
  
  # https://biocore-ntnu.github.io/pyranges/intersecting-ranges.html
  z_scores_filtered = z_scores.intersect(regions, how=how)

  if z_scores_filtered.df.empty: 
    return None, None, None, None, None

  z_scores_filtered.KBars = z_scores_filtered.KBars.replace({np.nan: None})

  return (
    z_scores_filtered.windowPositions.tolist(),
    z_scores_filtered.NBars.tolist(), 
    z_scores_filtered.NObserveds.tolist(), 
    z_scores_filtered.KBars.tolist(),
    z_scores_filtered.KObserveds.tolist()
  )

def compute_zscores_on_window(window, model, log=True):
  with pysam.TabixFile(model['mutations']) as mutations, pysam.FastaFile(model['genome']) as genome:
    N_bar, N_observed = compute_Nbar_Nobserved(window, model, mutations, genome, log) 
    K_bar, K_observed, M = compute_Kbar_Kobserved_M(window, model, mutations, genome)
  return N_bar, N_observed, K_bar, K_observed, M

def compute_expected_observed_counts(region, model, trustworthy_noncoding_regions_filename, window_stride, log=True):
  with pysam.TabixFile(model['mutations']) as mutations, pysam.FastaFile(model['genome']) as genome:
    windows = create_windows(model['windowSize'], window_stride, region, genome, region_contains_windows=True)
    N_bars, N_observeds = zip(*[
      compute_Nbar_Nobserved(window, model, mutations, genome, log) 
      for window in windows
    ])
    (
      SNV_positions_frequencies_CpG_positive, 
      SNV_positions_frequencies_CpG_negative
     ) = compute_SNV_positions_frequencies(mutations, genome, region, model)
    K_bars, K_observeds, Ms = zip(*[
      compute_Kbar_Kobserved_M(window, model, mutations, genome) 
      for window in windows
    ])

  # print([K_bar for K_bar in K_bars if not K_bar])
  # 1/0
  chromosome, start, end = unpack(region)

  (
    window_positions_trustworthy_noncoding_regions, 
    N_bars_trustworthy_noncoding_regions, 
    N_observeds_trustworthy_noncoding_regions,
    K_bars_trustworthy_noncoding_regions,
    K_observeds_trustworthy_noncoding_regions
  ) = filter_by_regions(
    windows,
    N_bars,
    N_observeds,
    K_bars,
    K_observeds, 
    regions=get_all_trustworthy_noncoding_regions(
      trustworthy_noncoding_regions_filename
    ), 
    how='containment'
  )
  (
    window_positions_exons, 
    N_bars_exons, 
    N_observeds_exons,
    K_bars_exons,
    K_observeds_exons
  ) = filter_by_regions(
    windows, 
    N_bars, 
    N_observeds, 
    K_bars, 
    K_observeds, 
    regions=get_canonical_exons(), 
    how=None
  )

  return {
    'region': region,
    'chromosome': chromosome,
    'start': start,
    'end': end,
    'windowPositions': [window['position'] for window in windows], 
    'windowRegions': [window['region'] for window in windows], 
    'windows': windows,
    'NBars': N_bars, 
    'NObserveds': N_observeds, 
    'snvCpGPositivePositions': pull_element(SNV_positions_frequencies_CpG_positive, index=0), 
    'snvCpGPositiveFrequencies': pull_element(SNV_positions_frequencies_CpG_positive, index=1), 
    'snvCpGNegativePositions': pull_element(SNV_positions_frequencies_CpG_negative, index=0),
    'snvCpGNegativeFrequencies': pull_element(SNV_positions_frequencies_CpG_negative, index=1),
    'KBars': K_bars, 
    'KObserveds': K_observeds, 
    'Ms': Ms,
    'windowPositionsTrustworthyNoncodingRegions': window_positions_trustworthy_noncoding_regions,
    'NBarsTrustworthyNoncodingRegions': N_bars_trustworthy_noncoding_regions,
    'NObservedsTrustworthyNoncodingRegions': N_observeds_trustworthy_noncoding_regions,
    'KBarsTrustworthyNoncodingRegions': K_bars_trustworthy_noncoding_regions,
    'KObservedsTrustworthyNoncodingRegions': K_observeds_trustworthy_noncoding_regions,
    'windowPositionsExons': window_positions_exons,
    'NBarsExons': N_bars_exons,
    'NObservedsExons': N_observeds_exons,
    'KBarsExons': K_bars_exons,
    'KObservedsExons': K_observeds_exons
  }

def parse_arguments(): 
  parser = argparse.ArgumentParser(description='')
  parser.add_argument('--region', type=str, help='')
  parser.add_argument('--model', type=str, help='')
  parser.add_argument('--window-stride', type=int, help='', dest='window_stride')
  parser.add_argument('--trustworthy-noncoding-regions', type=str, help='', dest='trustworthy_noncoding_regions')
  return parser.parse_args()

def test(): 
  args = parse_arguments()   

  from read_model import read_model
  model = read_model(args.model)

  print_json(compute_expected_observed_counts(
    args.region, 
    model, 
    args.trustworthy_noncoding_regions,
    args.window_stride
  ))

if __name__ == '__main__':
  test()