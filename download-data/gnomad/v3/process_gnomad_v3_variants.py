#!/usr/bin/env python3

from cyvcf2 import VCF
import argparse
import json
import itertools
import pandas as pd
import pyranges
import gzip

from pack_unpack import unpack, bed_to_sam_string

#%%
def parse_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--vep-annotation-file', type=str, dest='vep_annotation_file', help='')
    parser.add_argument('--filter-AC', default=0, type=int, dest='filter_AC', help='')
    parser.add_argument('--filter-AN', default=0, type=int, dest='filter_AN', help='')
    parser.add_argument('--filter-AF', default=0, type=int, dest='filter_AF', help='')    
    parser.add_argument('--gnomad-variant-file', type=str, dest='gnomad_variant_file', help='')
    #parser.add_argument('--coverage-file', type=str, dest='coverage_file', help='')
    parser.add_argument('--interval', type=str, dest='interval')
    parser.add_argument('--var-path', type=str, dest='var_path', help='')

    return parser.parse_args()
#%%

#%%
def read_variant_info(file):
    info_list = [] ## Initialize list 
    df = open(file, "r") ## Open the variant information file
    for line in df.readlines(): ## Iterate through each line of the information file
        info_list.append(line.rstrip()) ## Append variant info line into list 
    
    return info_list
    
#%%

#%%
def initialize_variant_dict(vep_annotations):
    ## Define list of key values for the variant's dictionary 
    variant_list = ['chrom', 'start', 'end', 'ref', 'alt', 'status', 'AC', 'AN', 'AF']
    variant_list.extend(vep_annotations)
    
    ## Use the vep annotation to define the keys for the dictionary
    variant_dict = dict.fromkeys(variant_list, '')    

    return variant_dict  

#%%

#%%
## Read in the variant VCF file with cyvcf2
## Iterate through each variant and get the necessary information
'''
1) variant coordinates
2) ref/alt alleles
3) vep annotations
'''
def get_variant_information(interval):
    
    ## Get all the arguments
    args = parse_arguments()
    
    # ## Read in coverage file 
    # coverage = pd.read_csv(args.coverage_file, sep="\t", header=0)
    # coverage.columns = "Chromosome Start End".split()
    
    ## Get vep_annotation names
    vep_annotations = read_variant_info(args.vep_annotation_file)
    
    ## Get the cyvcf2 object
    vcf = VCF(args.gnomad_variant_file)    
    
    ## Initialize list with full variant coordinate, allele, and vep information
    gnomad_v3_var = []
    
    #interval = "chr1:10100-10200"
    
    for variant in vcf(interval):
                
        ## Check to see if variants are present in the interval
        chrom = ''
        
        ## Get variant coordinates
        chrom, start, end = variant.CHROM, variant.start, variant.end
        # print(chrom, ":", start, "-", end, sep="")
        
        if len(chrom) == 0: 
            break
        
        ## Get the ref/alt allele
        ref, alt = variant.REF, variant.ALT
        
        ## Apply filters to determine if the variant should be kept or not
        keep_variant = filter_variants(variant, args.filter_AC, args.filter_AN, args.filter_AF)
        
        if keep_variant == True: 
        
            ## Get variant status
            status = "PASS"
            
            ## Get allele information
            ac = variant.INFO.get('AC')
            an = variant.INFO.get('AN')
            af = variant.INFO.get('AF')
            
            ## Get the vep annotations for the variant
            vep = variant.INFO.get('vep').split(';')
            
            ## Iterate through each isoform
            for isoform_num in range(0, len(vep)):
                ## Get the vep information for a specific isoform
                vep_isoform = vep[isoform_num].split('|')
            
                ## Initialize variant dictionary
                variant_dict = initialize_variant_dict(vep_annotations)
                
                ## Fill in dictionary with variant coordinates and ref/alt alleles
                variant_dict['chrom'], variant_dict['start'], variant_dict['end'] = chrom, start, end
                variant_dict['ref'], variant_dict['alt'] = list(itertools.chain(ref))[0], list(itertools.chain(alt))[0]
                variant_dict['status'], variant_dict['AC'], variant_dict['AN'], variant_dict['AF'] = status, ac, an, af
                
                ## Iterate through each isoform and add vep annotations to variant dictionary
                for annotation_num in range(0, len(vep_annotations)):

                    ## Define the vep annotation 
                    vep_annotation = vep_annotations[annotation_num]
                                    
                    ## Get the vep annotation value 
                    vep_value = vep_isoform[annotation_num]
                    
                    ## Fill in variant annotation value into dictionary
                    variant_dict[vep_annotation] = vep_value
                               
                ## Append to the final, full variant list 
                gnomad_v3_var.append(variant_dict)
    
    return gnomad_v3_var 
#%%

#%%
def filter_variants(variant, filter_AC, filter_AN, filter_AF):   

    ####################################################    
    ##### Perform filters on basic vcf info values #####
    ####################################################
    
    ## Remove variants whose filter value is not PASS
    ## https://brentp.github.io/cyvcf2/docstrings.html
    if str(variant.FILTER) != "None": 
        print('Sanity check failed: variant is not PASS... Removing...')
        return False
                    
    ## Filter variant based on alternate allele count 
    if filter_AC > 0: 
        if int(variant.INFO.get('AC')) < filter_AC: 
            print('Sanity check failed: variant allele acount (AC) below required threshold... Removing...')
            return False
    
    ## Filter variant based on alternate allele number 
    if filter_AN > 0: 
        if int(variant.INFO.get('AN')) < filter_AN: 
            print('Sanity check failed: variant allele number (AN) below required threshold... Removing...')
            return False
    
    ## Filter variant based on alterante allele fraction
    if filter_AF > 0: 
        if int(variant.INFO.get('AF')) < filter_AF:
            print('Sanity check failed: variant allele frequency (AF) below required threshold... Removing...')
            return False
        
    #################################
    ##### PERFORM SANITY CHECKS #####
    #################################
    
    ## Verify the variant is an SNV
    if variant.INFO.get('variant_type') == 'snv': 
        print('Sanity check failed: variant is not an SNV... Removing...')
        return False
    
    ## Verify reference and alternate allele is of length 1
    alt_allele = variant.ALT
    ref_allele = variant.REF
    
    if isinstance(alt_allele, list): 
        alt_allele = ''.join(map(str, alt_allele))
        
    if isinstance(ref_allele, list): 
        ref_allele = ''.join(map(str, ref_allele))
        
    alt_allele = list(itertools.chain(alt_allele))
    ref_allele = list(itertools.chain(ref_allele))
    
    if len(ref_allele) > 1: 
        print(variant.INFO.get('variant_type'))
        print(alt_allele)
        print(ref_allele)
        print('Sanity check failed: variant reference allele is of length > 1... Removing...')
        #return False
    
    if len(alt_allele) > 1: 
        print(variant.INFO.get('variant_type'))
        print(alt_allele)
        print(ref_allele)
        print('Sanity check failed: variant alternate allele is of length > 1... Removing...')
        #return False
    
    ## Verify ref/alt alleles are AGCT
    if ref_allele[0] not in ['A', 'C', 'T', 'G']: 
        print('Sanity check failed: ref allele not a valid nucleotide... Removing...')
        return False
    
    if alt_allele[0] not in ['A', 'C', 'T', 'G']: 
        print('Sanity check failed: alt allele not a valid nucleotide... Removing...')
        return False
    
    ## Verify that the ref and alt alleles are different
    if ref_allele[0] == alt_allele[0]: 
        print('Sanity check failed: ref and alt alleles are the same... Removing variant...')
        return False 
    
    ## Verify that the start/end coordinates of the SNV are 1 apart
    if (variant.end - variant.start) != 1: 
        print('Sanity check failed: variant start and end coordinates does not equal 1... Removing...')
        return False
    
    ############################################################
    ##### Remove variants within regions of lower coverage #####
    ############################################################
    ## Create pandas df from the variant coordinate
    # chrom, start, end = variant.CHROM, variant.start, variant.end
    # var = {'Chromosome': [chrom], 'Start': [start], 'End': [end]}
    # var = pd.DataFrame.from_dict(var)
    
    # ## Create pyranges object from the pandas df
    # cov, var = pyranges.PyRanges(coverage), pyranges.PyRanges(var)
    
    # ## Perform intersection
    # intersection = cov.intersect(var)
    # print(var)
    # print(intersection)
    
    # ## See if there is an intersection
    # if len(intersection) > 0: 
    #     print('Variant within insufficiently covered region.. Removing...')
    #     return False
    
    
    ## Print message saying the variant met all filtering criteria
    print("Variant met all filtering criteria...")
    return True

#%%

#%%
def process_gnomad_v3_variants():
    
    ## Read in the arguments
    args = parse_arguments()
    
    ## Perform checks on input arguemnts
    ## TODO
        
    annotated_vars = get_variant_information(args.interval)
        
    if annotated_vars: ## See if mutations were found within the interval
        if len(annotated_vars) > 0: ## See if any mutations met filtering criteria
            print(annotated_vars)
            
            ## Convert list of dictionaries to dataframe
            df = pd.DataFrame(annotated_vars)
            
            ## Select for necessary columns
            header = ['chrom', 'start', 'end', 'ref', 'alt', 'allele_count', 'gene_name', 'gene_id', 'consequence', 'feature_type', 'transcript_id', 'amino_acids']
            df = df[['chrom', 'start', 'end', 'ref', 'alt', 'AC', 'SYMBOL', 'Gene', 'Consequence', 'Feature_type', 'Feature', 'Amino_acids']]
            df.columns = header
        
            ## Convert list of dictionaries to json file
            variant_path = args.var_path \
                + '/intermediate_files/{}/'.format(str(df['chrom'].iloc[0])) \
                + 'gnomad_v3_variants_{}.bed'.format(args.interval)
            
            df.to_csv(variant_path, index=False, sep="\t", header=None)
            
            ## Write the header
            header = df.loc[[1]]
            
            header_path = args.var_path + '/header_{}'.format(str(df['chrom'].iloc[0]))
            header.to_csv(header_path, index=False, sep="\t")
            
            
    print("Job successfully ran for the interval: {}... Ready to merge...".format(args.interval))
    
#%%

if __name__ == '__main__':
    process_gnomad_v3_variants()
