import itertools 

bases = ['A', 'T', 'G', 'C']

def is_odd(kmer_size): 
  if kmer_size % 2 == 0: 
    raise ValueError('kmer size must be odd: {}'.format(kmer_size))

def compute_left_right(position, kmer_size, sequence_length): 
  is_odd(kmer_size)
  flank = int((kmer_size-1)/2)
  left = position - flank
  if left < 0: raise IndexError
  right = position + flank + 1
  if right > sequence_length: raise IndexError
  return left, right

def fetch_kmer_from_sequence(sequence, position, kmer_size): 
  left, right = compute_left_right(position, kmer_size, len(sequence))  
  return sequence[left:right]

def fetch_kmer_from_genome(genome, chromosome, position, kmer_size): 
  # provide the "reference" argument positionally to "get_reference_length": 
  # https://stackoverflow.com/a/24463222/6674256
  left, right = compute_left_right(position, kmer_size, genome.get_reference_length(chromosome))  
  # "fetch" API: https://pysam.readthedocs.io/en/latest/api.html?highlight=fasta#pysam.FastaFile
  return genome.fetch(chromosome, left, right).upper()

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

      

