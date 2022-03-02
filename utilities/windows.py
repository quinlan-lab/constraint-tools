import numpy as np 

# pysam API: 
# https://github.com/pysam-developers/pysam/blob/b82cbcae22c088e64fdb58f8acaf1e9773c7b088/pysam/libctabix.pyx
import pysam

from pack_unpack import unpack, pack
from colorize import print_json

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
  for window_center in np.arange(window_center_start, window_center_end+1, window_stride): 
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

def test(): 
  genome_filename = '/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/reference/grch38/hg38.analysisSet.fa.gz'
  region = 'chr1:15300-15310'
  with pysam.FastaFile(genome_filename) as genome: 
    print_json(create_windows(
      window_size=5, 
      window_stride=1, 
      region=region, 
      genome=genome, 
      region_contains_windows=True
    ))

if __name__ == '__main__': 
  test()