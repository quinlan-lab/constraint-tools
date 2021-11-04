# pysam API: 
# https://github.com/pysam-developers/pysam/blob/b82cbcae22c088e64fdb58f8acaf1e9773c7b088/pysam/libctabix.pyx
import pysam

from column_headings import fetch_column_heading_indices
from colorize import print_json
from kmer import fetch_kmer_from_genome, fetch_kmer_from_sequence, middle_base, middle_index, check_for_Ns

def fetch_SNVs(mutations, genome, region, meta):
  SNVs = []
  column_headings = [
    #'Chromosome'.lower(), 
    #'Start_Position'.lower(),
    #'End_Position'.lower(),
    #'Reference_Allele'.lower(), 
    #'Tumor_Seq_Allele1'.lower(), 
    #'Tumor_Seq_Allele2'.lower(), 
    #'Tumor_Sample_Barcode'.lower(),
    #'Variant_Type'.lower(), 
    #'ref_context' 
    'chrom', 'start', 'end', 'ref', 'alt', 'allele_count', 'variant_type'
  ]
  column_heading_indices = fetch_column_heading_indices(meta['mutations'])
  #print('column_heading_indices',column_heading_indices)
  for row in mutations.fetch(region=region, parser=pysam.asTuple()):
    mutation = {}
    for column_heading in column_headings: 
      try: 
        mutation[column_heading] = row[column_heading_indices[column_heading]]
      except KeyError: 
        continue
    #if mutation['reference_allele'] == '-': continue # exclude INS
    #print(mutation['start'])
    if mutation['ref'] == '-': continue # exclude INS
    #if len(mutation['reference_allele']) != 1: continue # exclude DEL, DNP, TNP, ONP
    if len(mutation['ref']) != 1: continue
    #if mutation['tumor_seq_allele2'] not in {'A', 'T', 'G', 'C'}: continue # exclude DEL
    if mutation['alt'] not in {'A', 'T', 'G', 'C'}: continue
    
    ###I don't think we need this block for this file format?
    #try: 
      #if mutation['reference_allele'] != mutation['tumor_seq_allele1']: # sanity check
    #    print_json(mutation)
    #    raise ValueError('Tumor allele 1 is not the same as reference allele') 
    #  del mutation['tumor_seq_allele1']
    #except KeyError:
    #  pass

    #start = int(mutation['start_position'])
    start = int(mutation['start'])
    #mutation['position'] = start - 1
    mutation['position'] = start
    #del mutation['start_position']
    del mutation['start']

    try: 
      #if start != int(mutation['end_position']): # sanity check for SNPs
      if start != int(mutation['end'])-1: # sanity check for SNPs
        print_json(mutation)
        #raise ValueError('SNPs must have equal start and end values') 
        raise ValueError('SNPs must start and end values 1 bp apart') 
      #del mutation['end_position']
      del mutation['end']
    except KeyError:
      pass

    #kmer = fetch_kmer_from_genome(genome, mutation['chromosome'], mutation['position'], meta['kmer_size'])
    kmer = fetch_kmer_from_genome(genome, mutation['chrom'], mutation['position'], meta['kmer_size'])
    mutation['kmer'] = kmer
    
    #if middle_base(kmer) != mutation['reference_allele'].upper(): # sanity check
    if middle_base(kmer) != mutation['ref'].upper(): # sanity check
      print_json(mutation)
      raise ValueError('Middle base of kmer does not match ref allele: {} {} {}'.format(
        middle_base(kmer), 
        kmer,
        #mutation['reference_allele']
        mutation['ref']
      )) 
    
    try:
      ref_context = mutation['ref_context'].upper() 
      kmer_from_maf = fetch_kmer_from_sequence(ref_context, middle_index(ref_context), meta['kmer_size'])        
      check_for_Ns(kmer_from_maf) # sanity check 
      if kmer_from_maf != kmer: # sanity check
        print_json(mutation)
        raise ValueError('kmer from maf does not match kmer inferred from fasta') 
      del mutation['ref_context']
    except KeyError:
      pass

    try: 
      #if mutation['variant_type'] != 'SNP': # sanity check
      if mutation['feature_type'] != 'SNP': # sanity check
        print_json(mutation)
        raise ValueError('variant_type is not SNP') 
      #del mutation['variant_type']
      del mutation['feature_type']
    except KeyError:
      pass

    #mutation['REF'] = mutation.pop('reference_allele')
    mutation['REF'] = mutation.pop('ref')
    #mutation['ALT'] = mutation.pop('tumor_seq_allele2')
    mutation['ALT'] = mutation.pop('alt')

    SNVs.append(mutation)
  return SNVs

