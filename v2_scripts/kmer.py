import itertools

from colorize import print_string_as_error 

bases = ['A', 'T', 'G', 'C']

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

def check_for_Ns(kmer): 
  if 'N' in kmer: 
    print_string_as_error('Please supply genomic sequences devoid of Ns')
    #raise ValueError
  return kmer

def fetch_kmer_from_sequence(sequence, position, kmer_size): 
  left, right = compute_left_right(position, kmer_size, len(sequence))  
  return check_for_Ns(sequence[left:right].upper())

def fetch_kmer_from_genome(genome, chromosome, position, kmer_size): 
  # provide the "reference" argument positionally to "get_reference_length": 
  # https://stackoverflow.com/a/24463222/6674256
  left, right = compute_left_right(position, kmer_size, genome.get_reference_length(chromosome))  
  # "fetch" API: https://pysam.readthedocs.io/en/latest/api.html?highlight=fasta#pysam.FastaFile
  return check_for_Ns(genome.fetch(chromosome, left, right).upper())

def middle_index(kmer): 
  kmer_size = len(kmer)
  is_odd(kmer_size)
  return int((kmer_size - 1)/2)
  
def middle_base(kmer): 
  return kmer[middle_index(kmer)]

def compute_kmers(kmer_size): 
  return [''.join(tup) for tup in itertools.product(bases, repeat=kmer_size)]

def alternate_bases(kmer): 
  return set(bases) - {middle_base(kmer)}
    
def CpG(kmer): 
  if len(kmer) == 1: return False 
  return (
      kmer[middle_index(kmer)] == 'C' and 
      kmer[middle_index(kmer)+1] == 'G'
  )

def not_CpG(kmer): 
  return not CpG(kmer)

def initialize_kmer_data(args): 
  return {
    kmer: {
      'CpG': CpG(kmer),
      'cohort_sequence_count': 0,
      'sequence_count': 0,
      'REF': middle_base(kmer),
      'ALT_counts': { alternate_base: 0 for alternate_base in alternate_bases(kmer) } 
    } for kmer in compute_kmers(args.kmer_size)
  }

def get_bases():
  return bases
