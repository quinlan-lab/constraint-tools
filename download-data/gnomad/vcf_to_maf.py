#!/usr/bin/env python3

#%%
import json
import parmap
import tqdm
import itertools
import argparse
import pysam
import pandas as pd

from kmer import fetch_kmer_from_genome
#%%

#%%
def get_maf_values(vcf_variant):
    
  
    ## Define variant column values to include from the vcf dataset
    maf_keys = ['chrom', 'start', 'end', 'ref', 'alt', 'VARIANT_CLASS']
    
    ## Select for relevant columns
    maf_variant = {maf_key: vcf_variant[maf_key] for maf_key in maf_keys}
    
    #################################
    ##### PERFORM SANITY CHECKS #####
    #################################
    
    ## Remove non-SNVs 
    if maf_variant['VARIANT_CLASS'] != 'SNV': 
        print('Variant is not an SNV... Removing variant...')
        return {}
    
    ## Verify that the reference and alternate allele is of length 1 (i.e. it is an SNV event)
    alt_allele = list(itertools.chain(maf_variant['alt']))
    ref_allele = list(itertools.chain(maf_variant['ref']))
    if len(ref_allele) > 1: 
        print("Sanity check failed: ref allele of length > 1 found... Removing variant...")
        return {}
    
    if len(alt_allele) > 1: 
        print("Sanity check failed: alt allele of length > 1 found... Removing variant...")
        return {}
    
    ## Verify that the ref/alt alleles are ACTG
    if ref_allele[0] not in ['A', 'C', 'T', 'G']: 
        print("Sanity check failed: ref allele (", ref_allele, ") not a valid nucleotide... Removing variant...", sep="")
        return {}
    
    if alt_allele[0] not in ['A', 'C', 'T', 'G']: 
        print("Sanity check failed: ref allele (", alt_allele, ") not a valid nucleotide... Removing variant...", sep="")
        return {}
        
    ## Check to see that the ref and alt alleles are different
    if ref_allele[0] == alt_allele[0]: 
        print("Sanity check failed: ref and alt alleles are the same... Removing variant...")
        return {}
    
    ## Check to see if the start/stop coordinates are 1 apart
    if (maf_variant['end'] - maf_variant['start']) != 1:
        print("Sanity check failed: start and stop coordinate separation does not equal 1... Removing variant...")
        return {}
    
    #####################################################################################    
    ##### Convert variant information values to those that are used in the MAF file #####
    #####################################################################################
    final_maf_variant = {}
    final_maf_variant['Chromosome'], final_maf_variant['Start_Position'], final_maf_variant['End_Position'] = \
        maf_variant['chrom'], maf_variant['start'], maf_variant['end']
    final_maf_variant['Reference_Allele'], final_maf_variant['Tumor_Seq_Allele1'], final_maf_variant['Tumor_Seq_Allele2'] = \
        ref_allele[0], ref_allele[0], alt_allele[0]
    final_maf_variant['Tumor_Sample_Barcode'] = 'germline'
    final_maf_variant['Variant_Type'] = 'SNP'
    
    # print(final_maf_variant)
    
    return final_maf_variant
#%%

#%%
## Get the kmer sequence for the variant #####
def get_kmer(final_maf_variant, genome, kmer_size):
        
    ## Get the reference kmer from the reference genome
    ref_kmer = fetch_kmer_from_genome(genome, final_maf_variant['Chromosome'], final_maf_variant['Start_Position'], kmer_size)
    
    ## Store the kmer value
    final_maf_variant['ref_context'] = ref_kmer
    
    return final_maf_variant
#%%

#%%
def parse_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--vcf', type=str, help='')
    parser.add_argument('--genome', type=str, help='')
    parser.add_argument('--kmer-size', type=int, dest='kmer_size', help='')
    parser.add_argument('--output', type=str, help='')
    
    return parser.parse_args()
#%%

#%%
def vcf_to_maf(vcf_variant):
    args = parse_arguments()
    
    if vcf_variant: 
        
        print(vcf_variant)
            
        ## Filter variants from VCF
        filtered_var = get_maf_values(vcf_variant[0])
        
        ## Get the reference kmer for the vcf variant \
        ## Read in the reference genome 
        genome = pysam.FastaFile(args.genome)
        
        ## Get the reference kmer and add it to the dictionary
        maf_var = get_kmer(filtered_var, genome, args.kmer_size) 
        #maf_var['Chromosome'] =  maf_var['Chromosome'].split('chr')[1]
        print(maf_var)
        
        return maf_var
#%%

#%%
## Read in the gnomad variant JSON file
file = "/scratch/ucgd/lustre-work/quinlan/u1240855/constraint-tools/data/gnomad/v3/gnomad_v3_variants.json"

vcf = json.load(open(file, "r"))
#%%

#%%
if __name__ == '__main__':
        
    args = parse_arguments()
    vcf = json.load(open(args.vcf, "r"))
    
    gnomad_maf = parmap.map(vcf_to_maf, vcf, pm_pbar=True)
                        
    ## Convert variant information to pandas dataframe
    df = pd.DataFrame(gnomad_maf)
    
    ## Sort pandas dataframe
    df.sort_values(by=['Chromosome', 'Start_Position'])

    ## Write gzip the sorted pandas dataset    
    df.to_csv(args.output, index=False, sep="\t")
    
#%%

#%%
