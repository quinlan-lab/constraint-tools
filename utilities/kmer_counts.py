from kmer import fetch_kmers
from snvs import fetch_SNVs, reduce_SNVs

from colorize import print_string_as_info

# using the "khmer" library might make this function a handful of times faster:
# https://khmer.readthedocs.io/en/v0.6.1-0/ktable.html
def compute_kmer_total_counts(region, genome, kmer_counts, args):
  for kmer in fetch_kmers(region, genome, args.kmer_size): 
    kmer_counts[kmer]['count'] += 1
  return kmer_counts

def compute_kmer_SNV_counts(region, mutations, genome, kmer_counts, args):
  print_string_as_info('Fetching SNVs in region {} and incrementing corresponding kmer counts...'.format(region))
  SNVs = fetch_SNVs(mutations, genome, region, args.__dict__)
  SNVs = reduce_SNVs(SNVs)
  for SNV in SNVs: 
    kmer_counts[SNV['kmer']]['ALTStateCounts'][SNV['ALTState']] += 1
  return kmer_counts

