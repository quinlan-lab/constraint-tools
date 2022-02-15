import itertools

bases = 'ACGT'

def is_odd(filter_size): 
  if filter_size % 2 == 0: 
    raise ValueError('filter size must be odd: {}'.format(filter_size))

def middle_index(kmer): 
  kmer_size = len(kmer)
  is_odd(kmer_size)
  return int((kmer_size - 1)/2)
  
def middle_base(kmer): 
  return kmer[middle_index(kmer)]

def get_bases():
  return bases

def get_alternate_bases(kmer): 
  return bases.replace(middle_base(kmer), '')
    
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

def contains_unspecified_bases(kmer): 
  return {'N', 'M', 'R'} & set(kmer)

def check_for_Ns(kmer): 
  # https://www.qmul.ac.uk/sbcs/iubmb/misc/naseq.html
  if contains_unspecified_bases(kmer):
    raise ValueError(f'kmer is: {kmer}\n') 
  return kmer

def fetch_kmer_from_sequence(sequence, position, kmer_size): 
  left, right = compute_left_right(position, kmer_size, len(sequence))  
  return sequence[left:right].upper()

def fetch_kmer_from_genome(genome, chromosome, position, kmer_size): 
  # provide the "reference" argument positionally to "get_reference_length": 
  # https://stackoverflow.com/a/24463222/6674256
  left, right = compute_left_right(position, kmer_size, genome.get_reference_length(chromosome))  
  # "fetch" API: https://pysam.readthedocs.io/en/latest/api.html?highlight=fasta#pysam.FastaFile
  return check_for_Ns(genome.fetch(chromosome, left, right).upper())

def compute_kmers(kmer_size): 
  return [''.join(tup) for tup in itertools.product(bases, repeat=kmer_size)]

def compute_possible_ALT_states_core(kmer, ALT_multiplicity):   
  return [
    '{' + ','.join(s) + '}' 
    for s in itertools.combinations(get_alternate_bases(kmer), ALT_multiplicity)
  ]

def compute_possible_ALT_states(kmer):   
  list_of_lists = [compute_possible_ALT_states_core(kmer, ALT_multiplicity) for ALT_multiplicity in [1, 2, 3]]
  return [item for sub_list in list_of_lists for item in sub_list]

def CpG(kmer): 
  if len(kmer) == 1: return False 
  return (
      kmer[middle_index(kmer)] == 'C' and 
      kmer[middle_index(kmer)+1] == 'G'
  )

def not_CpG(kmer): 
  return not CpG(kmer)

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

if __name__ == '__main__': 
  print(compute_kmers(3))
  print(compute_possible_ALT_states('AGCAT'))
