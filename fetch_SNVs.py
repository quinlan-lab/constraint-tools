# pysam API: 
# https://github.com/pysam-developers/pysam/blob/b82cbcae22c088e64fdb58f8acaf1e9773c7b088/pysam/libctabix.pyx
import pysam
from pyfaidx import Fasta

# TODO #1: for given maf, pull out the column numbers of the minimal column headers listed at: https://github.com/mskcc/vcf2maf/blob/main/data/minimalist_test_maf.tsv

# TODO #2: use this code to filter out indels, etc from maf
# https://github.com/petermchale/trfermikit/blob/207975c92b3207ca6d9e37955857a3f29329100b/filter-calls/find_SVs.py
# https://github.com/petermchale/trfermikit/blob/207975c92b3207ca6d9e37955857a3f29329100b/filter-calls/filterByUnitigSupport_annotate.py
# https://github.com/petermchale/trfermikit/blob/207975c92b3207ca6d9e37955857a3f29329100b/utilities/sv.py

def parse(region): 
  chromosome, start_end = region.split(':')
  start, end = map(lambda s: int(s.replace(',', '')), start_end.split('-'))
  return chromosome, start, end

def fetch_SNVs(maf, fasta, region): 
  with pysam.TabixFile(maf) as mutations, Fasta(fasta, as_raw=True) as genome: 

    print('genome: {}'.format(genome))
    print('sequence: {}'.format(genome['1'][57158:57196]))

    filtered_mutations = []
    for mutation in mutations.fetch(*parse(region), parser=pysam.asTuple()):
      chromosome, start, end, reference_allele, tumor_allele_1, tumor_allele_2 = mutation[1], mutation[2], mutation[3], mutation[7], mutation[8], mutation[9]
      kmer = genome[str(chromosome)][int(start)-1:int(end)].upper()
      if kmer != reference_allele: 
        print(kmer, reference_allele)
        1/0
      else: 
        print('pass')
      filtered_mutations.append(
        {
          'chromosome': chromosome,
          'start': start,
          'end': end,
          'reference_allele': reference_allele,
          'tumor_allele_1': tumor_allele_1,
          'tumor_allele_2': tumor_allele_2,
          # pass string NOT integer to __getitem__ method of Fasta object:
          'kmer': kmer
        }
      )
    return filtered_mutations

def test(): 
  mutations = fetch_SNVs(
    maf = "/scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools/data/icgc/mutations.sorted.maf.gz",
    fasta = "/scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools/data/reference/grch37/genome.fa.gz",
    region = '1:74,000-75,000'
  )
  for mutation in mutations: 
    print(mutation)

if __name__ == '__main__': 
  test()