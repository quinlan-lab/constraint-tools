# pysam API: 
# https://github.com/pysam-developers/pysam/blob/b82cbcae22c088e64fdb58f8acaf1e9773c7b088/pysam/libctabix.pyx
import pysam

from colorize import print_json, print_string_as_error
from kmer import fetch_kmer_from_genome, fetch_kmer_from_sequence, middle_base, middle_index, check_for_Ns

def fetch_SNVs(mutations, genome, region, meta):
  SNVs = []
  
  ## TO DO: Check that the mutations data has the appropriate fields 
    
  try: 
    fh = mutations.fetch(region=region, parser=pysam.asTuple())
  except ValueError: 
    print_string_as_error('No mutations found within {}...'.format(region))
    return 
  
  ## Iterate through each mutation found within the neutral region of interest
  for mutation in fh:
    
    print_string_as_error('got here')
    ## Initialize dictionary to store variant kmer data
    kmer_dict = {}
          
    ## Get the kmer sequence for the mutation
    kmer = fetch_kmer_from_genome(genome, mutation.chrom, mutation.start, meta['kmer_size'])
    kmer_dict['kmer'] = kmer
    
    if middle_base(kmer) != mutation.ref.upper(): # sanity check
      print_json(mutation)
      raise ValueError('Middle base of kmer does not match ref allele: {} {}'.format(
        middle_base(kmer), 
        mutation.ref.upper()
      )) 
    
    kmer_dict['REF'] = mutation.ref
    kmer_dict['ALT'] = mutation.alt
    
    SNVs.append(kmer_dict)
    
  return SNVs

