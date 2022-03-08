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

def fetch_subdict(d): 
  return {
    k: d[k]
    for k in ['AAA', 'AAC', 'AGA']
  }

class Args: 
  kmer_size = None
  
def fetch_args(): 
  args = Args()
  args.kmer_size = 3 
  return args

def test_aggregate_kmer_counts(): 
  from colorize import print_unbuffered, print_json
  print_unbuffered('')
  print_string_as_info('********* testing aggregate_kmer_counts function: ****************')
  print_unbuffered('')

  from kmer import initialize_kmer_counts_germline

  # https://realpython.com/introduction-to-python-generators/
  def fetch_counts(count_type, args, progress_bar):
    if count_type != 'kmerCounts': 
      raise ValueError

    x = initialize_kmer_counts_germline(args)
    x['AAA']['count'] += 1
    x['AAA']['count'] += 1
    x['AAA']['ALTStateCounts']['{C,T}'] += 1
    x['AAC']['count'] += 1
    x['AGA']['count'] += 1
    x['AGA']['ALTStateCounts']['{A}'] += 1
    print('x:')
    print_json(fetch_subdict(x))
    yield x 

    y = initialize_kmer_counts_germline(args)
    y['AAA']['count'] += 1
    y['AAA']['ALTStateCounts']['{C,T}'] += 1
    y['AGA']['count'] += 1
    y['AGA']['ALTStateCounts']['{C}'] += 1
    print('y:')
    print_json(fetch_subdict(y))
    yield y

    z = initialize_kmer_counts_germline(args)
    z['AGA']['count'] += 1
    z['AGA']['ALTStateCounts']['{C}'] += 1
    print('z:')
    print_json(fetch_subdict(z))
    yield z

  kmer_counts = aggregate_kmer_counts(fetch_counts, args=fetch_args(), progress_bar=None)
  print('aggregate x, y and z:') 
  print_json(fetch_subdict(kmer_counts))

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
  print('aggregate x, y and z:', singleton_counts)

if __name__ == '__main__': 
  test_aggregate_kmer_counts()
  test_aggregate_singleton_counts() 

