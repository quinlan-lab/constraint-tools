import functools

from kmer import add_kmer_counts_germline
from singleton import add_singleton_counts
from colorize import print_string_as_info

def add_counts(x, y, args): 
  d = {
    'kmerCounts': add_kmer_counts_germline(x['kmerCounts'], y['kmerCounts'], args),
  }
  if ('singletonCounts' in x.keys()) and ('singletonCounts' in y.keys()):
    d['singletonCounts'] = add_singleton_counts(x['singletonCounts'], y['singletonCounts'])
  return d

def aggregate_counts(fetch_counts, args, progress_bar): 
  return functools.reduce(
    lambda x, y: add_counts(x, y, args),
    fetch_counts(args, progress_bar)
  )

def fetch_subdict(d): 
  from singleton import make_serializable
  return {
    'kmerCounts': {
      k: d['kmerCounts'][k]
      for k in ['AAA', 'AAC', 'AGA']
    },
    'singletonCounts': make_serializable(d['singletonCounts'])
  }  

class Args: 
  kmer_size = None
  
def fetch_args(): 
  args = Args()
  args.kmer_size = 3 
  return args

def test_aggregate_counts(): 
  from colorize import print_unbuffered, print_json
  print_unbuffered('')
  print_string_as_info('********* testing aggregate_counts function: ****************')
  print_unbuffered('')

  from kmer import initialize_kmer_counts_germline
  from singleton import Singleton_Counts

  # https://realpython.com/introduction-to-python-generators/
  def fetch_counts(args, progress_bar):
    x = {
      'kmerCounts': initialize_kmer_counts_germline(args),
      'singletonCounts': Singleton_Counts()
    }
    x['kmerCounts']['AAA']['count'] += 1
    x['kmerCounts']['AAA']['count'] += 1
    x['kmerCounts']['AAA']['ALTStateCounts']['{C,T}'] += 1
    x['kmerCounts']['AAC']['count'] += 1
    x['kmerCounts']['AGA']['count'] += 1
    x['kmerCounts']['AGA']['ALTStateCounts']['{A}'] += 1
    x['singletonCounts'][5][3] += 1
    x['singletonCounts'][4][1] += 1
    x['singletonCounts'][5][1] += 1
    print('x:')
    print_json(fetch_subdict(x))
    yield x 

    y = {
      'kmerCounts': initialize_kmer_counts_germline(args),
      'singletonCounts': Singleton_Counts()
    }
    y['kmerCounts']['AAA']['count'] += 1
    y['kmerCounts']['AAA']['ALTStateCounts']['{C,T}'] += 1
    y['kmerCounts']['AGA']['count'] += 1
    y['kmerCounts']['AGA']['ALTStateCounts']['{C}'] += 1
    y['singletonCounts'][5][3] += 1
    y['singletonCounts'][6][2] += 1
    y['singletonCounts'][1][1] += 1
    print('y:')
    print_json(fetch_subdict(y))
    yield y

    z = {
      'kmerCounts': initialize_kmer_counts_germline(args),
      'singletonCounts': Singleton_Counts()
    }
    z['kmerCounts']['AGA']['count'] += 1
    z['kmerCounts']['AGA']['ALTStateCounts']['{C}'] += 1
    z['singletonCounts'][5][3] += 1 
    z['singletonCounts'][4][0] += 1
    z['singletonCounts'][4][1] += 1
    z['singletonCounts'][3][1] += 1
    print('z:')
    print_json(fetch_subdict(z))
    yield z

  counts = aggregate_counts(fetch_counts, args=fetch_args(), progress_bar=None)
  print('aggregate x, y and z:') 
  print_json(fetch_subdict(counts))

if __name__ == '__main__': 
  test_aggregate_counts()

