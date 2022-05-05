import argparse 
import collections
import numpy as np 

# pysam API: 
# https://github.com/pysam-developers/pysam/blob/b82cbcae22c088e64fdb58f8acaf1e9773c7b088/pysam/libctabix.pyx
import pysam

from get_p0s_p1s_p2s_p3s import get_p0s_p1s_p2s_p3s

def create_column_of_Ns_core(number_examples, p0_p1_p2_p3): 
  return np.random.choice(a=[0, 1, 2, 3], size=(number_examples, 1), p=p0_p1_p2_p3)

def fetch_distribution_N(window, model):
  number_examples = 100000
  create_column_of_Ns = lambda p0_p1_p2_p3: create_column_of_Ns_core(number_examples, p0_p1_p2_p3)
  with pysam.FastaFile(model['genome']) as genome:
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
def fetch_distribution_K(window, genome, model):
  # 1. compute observed SNV count for given window 
  # 2. pull out the singleton-count distribution for that SNV count (xs, ys)
  # 3. return xs, ys, and counts (where counts is the number of windows used to compute the distribution)
  pass

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

  from windows import create_windows 
  from colorize import print_json  

  with pysam.FastaFile(model['genome']) as genome:
    windows = create_windows(model['windowSize'], args.window_stride, args.region, genome, region_contains_windows=True)      
    window = windows[0]
  print_json(fetch_distribution_N(window, model))

if __name__ == '__main__':
  test()