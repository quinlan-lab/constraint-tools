#!/usr/bin/env python3

from cyvcf2 import VCF
import argparse
import json
import concurrent.futures
import numpy as np

#%%
def parse_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--intervals', type=str, help='')
    parser.add_argument('--vep_annotation_file', type=str, help='')
    parser.add_argument('--filter-AC', default=0, type=int, dest='filter_AC', help='')
    parser.add_argument('--filter-AN', default=0, type=int, dest='filter_AN', help='')
    parser.add_argument('--filter-AF', default=0, type=int, dest='filter_AF', help='')    
    parser.add_argument('--gnomad_variant_file', type=str, help='')
    parser.add_argument('--var_path', type=str, help='')

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
    
    ## Get vep_annotation names
    vep_annotations = read_variant_info(args.vep_annotation_file)
    
    ## Get the cyvcf2 object
    vcf = VCF(args.gnomad_variant_file)
    
    ## Define the interval of interest
    chromosome = 0
    start = 1
    end = 2
    
    chromosome = interval[chromosome]
    start = str(int(interval[start]))
    end = str(int(interval[end]))
    
    if chromosome == 23: 
        chromosome = "chrX"
    elif chromosome == 24: 
        chromosome = "chrY"
    else: 
        chromosome = 'chr' + str(int(chromosome))
        
    interval = chromosome + ":" + start + "-"  + end
    
    
    ## Initialize list with full variant coordinate, allele, and vep information
    gnomad_v3_var = []
    
    for variant in vcf(interval): 
        ## Check to see if variants are present in the interval
        chrom = ''
        
        ## Get variant coordinates
        chrom, start, end = variant.CHROM, variant.start, variant.end
        #print(chrom, ":", start, "-", end, sep="")
        
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
                variant_dict['chrom'], variant_dict['start'], variant_dict['end'], variant_dict['ref'], variant_dict['alt'] = chrom, start, end, ref, alt
                variant_dict['status'], variant_dict['AC'], variant_dict['AN'], variant_dict['AF'] = status, ac, an, af
                
                ## Iterate through each isoform and add vep annotations to variant dictionary
                for annotation_num in range(0, len(vep_annotations)):
                    ## Define the vep annotation 
                    vep_annotation = vep_annotations[annotation_num]
                    #print("Getting the variant's VEP annotation: ", vep_annotation, "...", sep="")
                                    
                    ## Get the vep annotation value 
                    vep_value = vep_isoform[annotation_num]
                    
                    ## Fill in variant annotation value into dictionary
                    variant_dict[vep_annotation] = vep_value
                               
                ## Append to the final, full variant list 
                gnomad_v3_var.append([variant_dict])
    
        return gnomad_v3_var    
#%%

#%%
def filter_variants(variant, filter_AC, filter_AN, filter_AF):   
    
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
        
    if filter_AF > 0: 
        if int(variant.INFO.get('AF')) < filter_AF:
            print('Sanity check failed: variant allele frequency (AF) below required threshold... Removing...')
            return False
    
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
        
    ## Get the intervals to split gnomad variants
    intervals = np.loadtxt(args.intervals, delimiter="\t")
    
    ## Execute multiprocessing over all the intervals
    with concurrent.futures.ProcessPoolExecutor() as executor: 
        annotated_vars = executor.map(get_variant_information, intervals)
    
    ## Combine data
    my_list = []
    
    for var in annotated_vars: 
        if var:  ## Remove null values due to the lack of variants within soecfied intervals 
            my_list.append(var)

    ## Convert list of dictionaries to json file
    variant_path = args.var_path + '/gnomad_v3_variants.json'
    with open(variant_path, 'w') as fh:
        json.dump(my_list, fh, indent=2)
    
#%%

if __name__ == '__main__':
    process_gnomad_v3_variants()
