from ast import Assert
import numpy as np 

# pysam API: 
# https://github.com/pysam-developers/pysam/blob/b82cbcae22c088e64fdb58f8acaf1e9773c7b088/pysam/libctabix.pyx
import pysam

from pack_unpack import unpack, pack
from colorize import print_json, print_string_as_error

def is_odd(filter_size): 
  if filter_size % 2 == 0: 
    raise ValueError('filter size must be odd: {}'.format(filter_size))

def compute_left_right(position, filter_size, sequence_length, offset='zero_offset'): 
  is_odd(filter_size)
  flank = int((filter_size-1)/2)
  left = position - flank
  if left < 0: raise IndexError
  if offset == 'zero_offset':
    right = position + flank + 1
  elif offset == 'unit_offset':
    right = position + flank
  else: 
    raise ValueError
  if right > sequence_length: raise IndexError
  return left, right

# if the observed "time series" comes from overlapping windows (window_stride < window_size),
# then the time series is auto-correlated, 
# necessitating the use of a HMM to model the time series
def create_windows(window_size, window_stride, region, genome, region_contains_windows): 
  windows = []
  chromosome, region_start, region_end = unpack(region)
  window_center_start = region_start + window_size//2 if region_contains_windows else region_start
  window_center_end = region_end - window_size//2 if region_contains_windows else region_end
  for window_center in np.arange(window_center_start, window_center_end, window_stride): 
    # provide the "reference" argument positionally to "get_reference_length": 
    # https://stackoverflow.com/a/24463222/6674256
    window_start, window_end = compute_left_right(
      window_center, 
      window_size, 
      genome.get_reference_length(chromosome), 
      offset='unit_offset'
    )  
    windows.append({
      'region': pack(chromosome, window_start, window_end),
      'position': int(window_center) # https://stackoverflow.com/a/50916741/6674256
    })
  return windows

def create_window(region): 
  chromosome, region_start, region_end = unpack(region)
  window_size = region_end - region_start
  window_center = region_start + window_size//2 
  return {
      'region': pack(chromosome, region_start, region_end),
      'position': int(window_center) # https://stackoverflow.com/a/50916741/6674256
  }

def test_create_windows(): 
  from colorize import print_string_as_info, print_string_as_info_dim 
  from kmer import truncate, fetch_kmer_from_genome
  genome_filename = '/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/reference/grch38/hg38.analysisSet.fa.gz'
  region = 'chr1:15300-15310'
  kmer_size = window_size = 5
  with pysam.FastaFile(genome_filename) as genome: 
    sequence = genome.fetch(*unpack(region))    
    print_string_as_info(f"Sequence for region {region}:")
    print_string_as_info_dim(truncate(sequence))
    windows = create_windows(
      window_size=window_size, 
      window_stride=2,
      region=region, 
      genome=genome, 
      region_contains_windows=True
    )
    print_json(windows)
    expected_kmers = ['GGCAG', 'CAGCT', 'GCTTG']
    for i, window in enumerate(windows): 
      chromosome, _, _ = unpack(region)
      kmer = fetch_kmer_from_genome(genome, chromosome, window['position'], kmer_size)
      print(kmer)
      try: 
        assert kmer == expected_kmers[i]
      except AssertionError: 
        print_string_as_error('kmer is not expected!')
    try: 
      assert (i + 1) == len(expected_kmers)
    except AssertionError: 
      print_string_as_error('number of kmers is not the expected number!')
    

if __name__ == '__main__': 
  test_create_windows()