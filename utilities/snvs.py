from collections import defaultdict 
import functools
import color_traceback 

# pysam API: 
# https://github.com/pysam-developers/pysam/blob/b82cbcae22c088e64fdb58f8acaf1e9773c7b088/pysam/libctabix.pyx
import pysam

from column_headings import fetch_column_headings_indices
from colorize import print_json
from kmer import fetch_kmer_from_genome, middle_base

def fetch_SNVs(mutations, genome, region, meta, number_chromosomes_min=0):
  SNVs = []
  column_headings, heading_to_index = fetch_column_headings_indices(meta['mutations'])
  for row in mutations.fetch(region=region, parser=pysam.asTuple()):
    row_dict = {}
    for column_heading in column_headings: 
      try: 
        row_dict[column_heading] = row[heading_to_index[column_heading]]
      except KeyError: 
        continue

    if len(row_dict['REF']) > 1: continue # exclude DEL, DNP, TNP, ONP
    if len(row_dict['ALT']) > 1: continue # exclude INS 
    
    if int(row_dict['start']) != int(row_dict['end']) - 1: # sanity check for SNPs
      print_json(mutation)
      raise ValueError('SNPs must start and end 1 bp apart') 

    if row_dict['variant_type'] not in {'snv', 'multi-snv', 'mixed'}: # sanity check for SNPs
      print_json(mutation)
      raise ValueError("SNPs must be one of 'snv', 'multi-snv' or 'mixed'") 

    # filter out SNVs that are genotyped in too few chromosomes
    if int(row_dict['number_chromosomes']) < number_chromosomes_min: continue 

    mutation = {
      'chromosome': row_dict['chromosome'],
      'position': int(row_dict['start']),
      'REF': row_dict['REF'],
      'ALT': row_dict['ALT'],
      'number_ALT_chromosomes': int(row_dict['number_ALT_chromosomes']),
      'number_chromosomes': int(row_dict['number_chromosomes']),
    }

    kmer = fetch_kmer_from_genome(genome, mutation['chromosome'], mutation['position'], meta['kmer_size'])
    if middle_base(kmer) != mutation['REF'].upper(): # sanity check
      print_json(mutation)
      raise ValueError(f"Middle base of kmer does not match ref allele: {kmer} {mutation['REF']}") 
    mutation['kmer'] = kmer 

    SNVs.append(mutation)
  return SNVs

def combine_SNVs(x, y):
  # sanity check: 
  if x['kmer'] != y['kmer']:
    print_json(x)
    print_json(y)
    raise ValueError("SNVs in the same (position) group should have identical kmers!")  

  return {
    'kmer': x['kmer'],
    'ALT': x['ALT'] + y['ALT']
  }

def create_ALT_state(SNV_reduced): 
  sorted_alleles = sorted(SNV_reduced['ALT'])
  return { 
    'kmer': SNV_reduced['kmer'],
    'ALTState': '{' + ','.join(sorted_alleles) + '}'
  }

def reduce_SNVs(SNVs): 
  # group SNVs by position:
  positions_to_SNVs = defaultdict(lambda: []) 
  for SNV in SNVs: 
    positions_to_SNVs[SNV['position']].append({
      'kmer': SNV['kmer'],
      'ALT': [SNV['ALT']]
    })

  # reduce all SNVs in a group to a single SNV:
  SNVs_reduced = []
  for SNVs_at_position in positions_to_SNVs.values(): 
    SNVs_reduced.append(functools.reduce(combine_SNVs, SNVs_at_position))

  return [create_ALT_state(SNV_reduced) for SNV_reduced in SNVs_reduced]

def print_variants(vcf_filename, region): 
  from cyvcf2 import VCF
  vcf = VCF(vcf_filename)
  for variant in vcf(region):
    print_json({
      'chromosome': variant.CHROM,
      'start': variant.start,
      'end': variant.end,
      # Possible values are "snv", "indel", "multi-snv", "multi-indel", or "mixed" :
      'variant_type': variant.INFO.get('variant_type'),
      'REF': variant.REF,
      'ALT': variant.ALT[0],
      # Total number of alternate alleles observed at variant's locus :
      'number_ALT': variant.INFO.get('n_alt_alleles'),
      ##INFO=<ID=AC,Number=A,Type=Integer,Description="Alternate allele count">
      'number_ALT_chromosomes': variant.INFO.get('AC'),
      ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles">
      'number_chromosomes': variant.INFO.get('AN'),
    })

def test_fetch_SNVs(): 
  mutations_filename = '/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/gnomad/v3/variants/gnomad_v3.sorted.tsv.gz'
  genome_filename = '/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/reference/grch38/hg38.analysisSet.fa.gz'
  neutral_region = 'chr1:15363-15768'
  meta = { 
    'mutations': mutations_filename, 
    'kmer_size': 3
  }
  number_chromosomes_min = 150000
  with pysam.TabixFile(mutations_filename) as mutations, pysam.FastaFile(genome_filename) as genome: 
    SNVs = fetch_SNVs(mutations, genome, neutral_region, meta)
    SNVs_filtered = fetch_SNVs(mutations, genome, neutral_region, meta, number_chromosomes_min)

    print('SNVs at positions where more than one ALT allele is found in the cohort:')
    print_json([SNV for SNV in SNVs if (SNV['position'] == 15556) or (SNV['position'] == 15658)])

    print(f'number of SNVs: {len(SNVs)}')    
    print(f'number of SNVs genotyped in more than {number_chromosomes_min} chromosomes: {len(SNVs_filtered)}')

  vcf_filename = '/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/gnomad/v3/variants/gnomad_v3_chr1.vcf.gz'

  print('')
  print('query vcf directly (instead of tsv, which is a processed version of vcf)...')
  print('two adjacent polymorphic sites, the second of which has the following ALT allele profiles... ')
  print('https://gnomad.broadinstitute.org/variant/1-15557-G-A?dataset=gnomad_r3')
  print('https://gnomad.broadinstitute.org/variant/1-15557-G-C?dataset=gnomad_r3')
  print('https://gnomad.broadinstitute.org/variant/1-15557-G-T?dataset=gnomad_r3')
  print_variants(vcf_filename, 'chr1:15556-15557')

  print('query vcf directly (instead of tsv, which is a processed version of vcf) at another location with multiple ALTs...')
  print('https://gnomad.broadinstitute.org/variant/1-15659-G-C?dataset=gnomad_r3')
  print('https://gnomad.broadinstitute.org/variant/1-15659-G-A?dataset=gnomad_r3')
  print_variants(vcf_filename, 'chr1:15658-15659')

  return SNVs

def test_reduce_SNVs(SNVs): 
  print('Sites where more than one ALT allele is segregating: ')
  print_json([SNV for SNV in reduce_SNVs(SNVs) if len(SNV['ALTState']) > 3])

def test(): 
  SNVs = test_fetch_SNVs()
  test_reduce_SNVs(SNVs)

if __name__ == '__main__': 
  test()
