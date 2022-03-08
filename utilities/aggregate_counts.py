import functools

from kmer import combine_kmer_counts_germline
from singleton import combine_singleton_counts
from colorize import print_string_as_info

def aggregate_kmer_counts(fetch_counts, args, progress_bar): 
  print_string_as_info('Combining kmer counts...')
  return functools.reduce(
    lambda x, y: combine_kmer_counts_germline(x, y, args), 
    fetch_counts('kmerCounts', args, progress_bar)
  )

def aggregate_singleton_counts(fetch_counts, args, progress_bar): 
  print_string_as_info('Combining singleton counts...')
  return functools.reduce(
    combine_singleton_counts, 
    fetch_counts('singletonCounts', args, progress_bar)
  )

def aggregate_counts(fetch_counts, args, kmerCounts_log, singletonCounts_log): 
  kmer_counts = aggregate_kmer_counts(fetch_counts, args, kmerCounts_log)
  singleton_counts = aggregate_singleton_counts(fetch_counts, args, singletonCounts_log)
  return kmer_counts, singleton_counts

def test_aggregate_singleton_counts(): 
  from colorize import print_unbuffered
  print_unbuffered('')
  print_string_as_info('********* testing aggregate_singleton_counts function: ****************')
  print_unbuffered('')

  from singleton import Singleton_Counts

  # https://realpython.com/introduction-to-python-generators/
  def fetch_counts(count_type, args, progress_bar):
    if count_type != 'singletonCounts': 
      raise ValueError

    x = Singleton_Counts()
    x[5][3] += 1
    x[4][1] += 1
    x[5][1] += 1
    print('x:', x)
    yield x 

    y = Singleton_Counts()
    y[5][3] += 1
    y[6][2] += 1
    y[1][1] += 1
    print('y:', y)
    yield y 

    z = Singleton_Counts()
    z[5][3] += 1 
    z[4][0] += 1
    z[4][1] += 1
    z[3][1] += 1
    print('z:', z)
    yield z 

  singleton_counts = aggregate_singleton_counts(fetch_counts, args=None, progress_bar=None)
  print('sum of x, y and z:', singleton_counts)

if __name__ == '__main__': 
  test_aggregate_singleton_counts() 

