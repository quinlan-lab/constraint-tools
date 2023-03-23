import itertools
import color_traceback 

from colorize import (
  print_string_as_error,
  print_string_as_info,
  print_string_as_info_dim
)
from windows import (
  is_odd, 
  compute_left_right,
  create_windows
)
from pack_unpack import unpack

bases = 'ACGT'
complement = {
  'A': 'T', 
  'C': 'G', 
  'G': 'C', 
  'T': 'A'
}

def get_reverse_complement(kmer): 
  return ''.join(complement[base] for base in reversed(kmer))

def get_complement(ALT_state): 
  ALT_alleles = ALT_state.strip('{}').split(',')
  ALT_alleles_complemented = (complement[ALT_allele] for ALT_allele in ALT_alleles)
  return '{' + ','.join(sorted(ALT_alleles_complemented)) + '}'

def check_odd(kmer_size): 
  if not is_odd(kmer_size):     
    raise ValueError('kmer size must be odd: {}'.format(kmer_size))

def middle_index(kmer): 
  kmer_size = len(kmer)
  check_odd(kmer_size)
  return int((kmer_size - 1)/2)
  
def middle_base(kmer): 
  return kmer[middle_index(kmer)]

def get_alternate_bases(kmer): 
  return bases.replace(middle_base(kmer), '')
    
def contains_unspecified_bases(kmer): 
  # https://www.qmul.ac.uk/sbcs/iubmb/misc/naseq.html
  return {'N', 'M', 'R'} & set(kmer)

def fetch_kmer_from_genome_core(genome, chromosome, position, kmer_size): 
  # provide the "reference" argument positionally to "get_reference_length": 
  # https://stackoverflow.com/a/24463222/6674256
  left, right = compute_left_right(position, kmer_size, genome.get_reference_length(chromosome))  
  # "fetch" API: https://pysam.readthedocs.io/en/latest/api.html?highlight=fasta#pysam.FastaFile
  return genome.fetch(chromosome, left, right).upper()

def fetch_kmer_from_genome(genome, chromosome, position, kmer_size): 
  check_odd(kmer_size)
  return fetch_kmer_from_genome_core(genome, chromosome, position, kmer_size)

def truncate(string, length=100):   
  return string[:length] + '...' if len(string) > length else string

def fetch_kmers(region, genome, kmer_size, log=True): 
  # "fetch" API: https://pysam.readthedocs.io/en/latest/api.html?highlight=fasta#pysam.FastaFile
  # Note that fetch(region=region) does not work if the coordinates in "region" contains commas
  # Workaround is to parse "region" into "chromosome", "start", "end": 
  sequence = genome.fetch(*unpack(region))    
  if log: 
    print_string_as_info(f"Sequence for region {region}:")
    print_string_as_info_dim(truncate(sequence))

  windows = create_windows(
    window_size = kmer_size, 
    window_stride = 1, 
    region = region, 
    genome = genome, 
    # c.f. fetch_SNVs
    # region_contains_windows = True
    region_contains_windows = False
  )
  number_of_sites = 0 
  number_of_sites_containing_unspecified_bases = 0
  if log: 
    print_string_as_info('Iterating over region {} ...'.format(region))
  for window in windows: 
    number_of_sites += 1
    chromosome, _, _ = unpack(region)
    kmer = fetch_kmer_from_genome(genome, chromosome, window['position'], kmer_size)
    if contains_unspecified_bases(middle_base(kmer)): 
      number_of_sites_containing_unspecified_bases += 1
    if contains_unspecified_bases(kmer): continue 
    yield kmer 
    
  if log: 
    print_string_as_info_dim(f'Interrogated {number_of_sites}/{len(sequence)} sites in region {region}')
    print_string_as_info_dim(
      f'Number of sites in {region} containing unspecified bases: '
      f'{number_of_sites_containing_unspecified_bases}/{len(sequence)}'
    )

def compute_kmers(kmer_size): 
  return [''.join(tup) for tup in itertools.product(bases, repeat=kmer_size)]

def compute_possible_ALT_states_core(kmer, ALT_multiplicity):   
  return [
    '{' + ','.join(sorted(s)) + '}' 
    for s in itertools.combinations(get_alternate_bases(kmer), ALT_multiplicity)
  ]

def compute_possible_ALT_states(kmer):   
  list_of_lists = [compute_possible_ALT_states_core(kmer, ALT_multiplicity) for ALT_multiplicity in [1, 2, 3]]
  return [item for sub_list in list_of_lists for item in sub_list]

def C_followed_by_G(kmer): 
  return (
      kmer[middle_index(kmer)] == 'C' and 
      kmer[middle_index(kmer)+1] == 'G'
  )

def A_followed_by_T(kmer): 
  return (
      kmer[middle_index(kmer)] == 'A' and 
      kmer[middle_index(kmer)+1] == 'T'
  )

def CpG(kmer):
  if len(kmer) == 1: return False 
  return C_followed_by_G(kmer) or C_followed_by_G(get_reverse_complement(kmer))

def ApT(kmer):
  if len(kmer) == 1: return False 
  return A_followed_by_T(kmer) or A_followed_by_T(get_reverse_complement(kmer))

def not_CpG(kmer): 
  return not CpG(kmer)

def not_ApT(kmer): 
  return not ApT(kmer)

def initialize_kmer_counts_cancer(args): 
  return {
    kmer: {
      'CpG': CpG(kmer),
      'cohort_sequence_count': 0,
      'sequence_count': 0,
      'REF': middle_base(kmer),
      'ALT_counts': { alternate_base: 0 for alternate_base in get_alternate_bases(kmer) } 
    } for kmer in compute_kmers(args.kmer_size)
  }

def initialize_kmer_counts_germline(args): 
  return {
    kmer: {
      'CpG': CpG(kmer),
      'count': 0,
      'REF': middle_base(kmer),
      'ALTStateCounts': { ALT_state: 0 for ALT_state in compute_possible_ALT_states(kmer) },
    } for kmer in compute_kmers(args.kmer_size)
  }

def add_kmer_counts_germline(x, y, args):
  return {
    kmer: {
      'CpG': CpG(kmer),
      'count': int(x[kmer]['count']) + int(y[kmer]['count']),
      'REF': middle_base(kmer),
      'ALTStateCounts': {
        ALT_state: int(x[kmer]['ALTStateCounts'][ALT_state]) + int(y[kmer]['ALTStateCounts'][ALT_state])
        for ALT_state in compute_possible_ALT_states(kmer) 
      },
    } for kmer in compute_kmers(args.kmer_size)
  }

def test_get_reverse_complement(kmer): 
  print(f"reverse complement of {kmer} is {get_reverse_complement(kmer)}")

def test_get_complement(kmer):
  print('')
  print_string_as_info('**** testing get_complement(ALT_state) ***** ')
  for ALT_state in compute_possible_ALT_states(kmer): 
    print(kmer, ALT_state, get_complement(ALT_state))

def test_fetch_kmers(): 
  print('')
  print_string_as_info('**** testing fetch_kmers(...) ***** ')

  import pysam 

  genome_filename = '/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/reference/grch38/hg38.analysisSet.fa.gz'
  region = 'chr1:15300-15310'
  region_with_flanks = 'chr1:15298-15312'
  kmer_size = 5
  expected_kmers = ['GAGGC', 'AGGCA', 'GGCAG', 'GCAGC', 'CAGCT', 'AGCTT', 'GCTTG', 'CTTGC', 'TTGCC', 'TGCCT'] 
  with pysam.FastaFile(genome_filename) as genome: 
    sequence_with_flanks = genome.fetch(*unpack(region_with_flanks))    
    print_string_as_info(f'sequence for region {region_with_flanks}')
    print_string_as_info_dim(sequence_with_flanks)

    number_fetched_kmers = 0
    for i, kmer in enumerate(fetch_kmers(region, genome, kmer_size=kmer_size)): 
      number_fetched_kmers += 1
      print_string_as_info_dim(f'kmer: {kmer}')
      try: 
        assert kmer == expected_kmers[i]
      except AssertionError: 
        print_string_as_error('fetched kmer is not expected!')
    try: 
      assert number_fetched_kmers == len(expected_kmers)
    except AssertionError: 
      print_string_as_error('number fetched kmers differs from expected number of kmers!')
  print_string_as_info_dim('Tests passed')

def test_fetch_kmer_from_genome(): 
  print('')
  print_string_as_info('**** testing fetch_kmer_from_genome(...) ***** ')

  import pysam 

  genome_filename = '/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/reference/grch38/hg38.analysisSet.fa.gz'
  positions = [63442000, 63442001, 63442002, 63442003, 63442004]
  expected_kmers = ['A', 'G', 'A', 'G', 'T']
  observed_kmers = []
  with pysam.FastaFile(genome_filename) as genome: 
    for position in positions:
      kmer = fetch_kmer_from_genome(genome, chromosome='chr20', position=position, kmer_size=1) 
      print(f'nucleotide at position {position} is {kmer}')
      observed_kmers.append(kmer)
    try: 
      assert observed_kmers == expected_kmers
    except AssertionError:
      print_string_as_error('nucleotides are unexpected!')
  print_string_as_info_dim('Tests passed')

if __name__ == '__main__': 
  # print(compute_kmers(3))
  print("possible ALT states of 'AGTAT':", compute_possible_ALT_states('AGTAT'))
  print("possible ALT states of 'AGTAT' of multiplicity 2:", compute_possible_ALT_states_core(kmer='AGTAT', ALT_multiplicity=2))
  print("CpG('AACGT'):", CpG('AACGT'))
  print("CpG('ACGTT'):", CpG('ACGTT'))
  print("CpG('AACCT'):", CpG('AACCT'))
  print("CpG('AGGTT'):", CpG('AGGTT'))
  print("CpG('CGT'):", CpG('CGT'))
  print("CpG('AGT'):", CpG('AGT'))
  test_get_reverse_complement('AGCGT')
  test_get_reverse_complement('ACGCT')
  test_get_complement('AGTAT')
  test_get_complement('AGCGT')
  test_fetch_kmers()
  test_fetch_kmer_from_genome()
  