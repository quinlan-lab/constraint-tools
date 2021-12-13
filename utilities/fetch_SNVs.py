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
  print("This is the list of mutations for the region: ")  
  for mutation in fh:
    print(mutation)

    ## Initialize dictionary to store variant kmer data
    kmer_dict = {}
    
    ## Humanize varialbles (TODO: make this more robust)
    chrom = 0
    start = 1
    ref = 3
    alt = 4
    ac = 5
          
    ## Get the kmer sequence for the mutation
    kmer = fetch_kmer_from_genome(genome, mutation[chrom], int(mutation[start]), meta['kmer_size'])
    kmer_dict['kmer'] = kmer
    
    ## Sanity check: see if the middle base of the kmer does not equal the reference genome
    if middle_base(kmer) != mutation[ref].upper(): # sanity check
      print_json(mutation)
      raise ValueError('Middle base of kmer does not match ref allele: {} {}'.format(
        middle_base(kmer), 
        mutation[ref].upper()
      )) 
    
    ## Get the reference and alterante allele
    kmer_dict['REF'] = mutation[ref]
    kmer_dict['ALT'] = mutation[alt]
    kmer_dict['allele_count'] = mutation[ac]
    
    SNVs.append(kmer_dict)
    
  return SNVs

