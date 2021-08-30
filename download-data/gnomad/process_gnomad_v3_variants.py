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
    # parser.add_argument('--inclusion_variants_file', default=[], type=str, help='') ## Fix these defaults
    # parser.add_argument('--exclusion_variants_file', default=[], type=str, help='') ## Fix these defaults
    # parser.add_argument('--filter_variant_consequence', action=argparse.BooleanOptionalAction, help='If true, must provide inclusion and exclusion variant file...')
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
    
    # ## Define inclusion/exclusion variants 
    # ## Refer to https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html for variant consequences ordered by severity
    # inclusion_variants = read_variant_info(args.inclusion_variants_file)
    # exclusion_variants = read_variant_info(args.exclusion_variants_file)
    
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
        ## Get variant coordinates
        chrom, start, end = variant.CHROM, variant.start, variant.end
        #print(chrom, ":", start, "-", end, sep="")
        
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
                    
                    # if vep_annotation == "VARIANT_CLASS": 
                    #     if vep_isoform != "SNV": 
                    #             print('Sanity check failed: variant is not classified as an SNV... Removing...')      
                    #             continue
                    
                    ## Fill in variant annotation value into dictionary
                    variant_dict[vep_annotation] = vep_value
                
                
                ## Append to the final, full variant list 
                gnomad_v3_var.append([variant_dict])
    
        return gnomad_v3_var    
#%%

#%%
def filter_variants(variant, filter_AC, filter_AN, filter_AF):   
    
    # ## Remove variants in the exclusion list 
    # if (consequence in exclusion_variants) and (filter_variant_consequence == True): 
    #     print("Removing ", consequence, " variant...", sep="")
    #     return False
        
    # ## Include variants in the inclusion list
    # if (consequence in inclusion_variants) and (filter_variant_consequence == True): 
    #     print("Including ", consequence, " variant", sep="")
    #     print("Checking to see if the variant meets additional filtering criteria...")
        
    # ## Remove variants not found in either list 
    # if (consequence not in inclusion_variants and consequence not in exclusion_variants) and (filter_variant_consequence == True):
    #     print("Variant of consequence: ", consequence, " not in the inclusion or exclusion list... Removing...", sep="")
    #     return False
    
    ## Remove variants whose filter value is not PASS
    ## https://brentp.github.io/cyvcf2/docstrings.html
    if variant.FILTER != "None": 
        print('Sanity check failed: variant is not PASS... Removing...')
        return False
                    
    ## Filter variant based on alternate allele count 
    if filter_AC > 0: 
        if variant.INFO.get('AC') < filter_AC: 
            print('Sanity check failed: variant allele acount (AC) below required threshold... Removing...')
            return False
    
    ## Filter variant based on alternate allele number 
    if filter_AN > 0: 
        if variant.INFO.get('AN') < filter_AN: 
            print('Sanity check failed: variant allele number (AN) below required threshold... Removing...')
            return False
        
    if filter_AF > 0: 
        if variant.INFO.get('AF') < filter_AF:
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

    ## Process vep annotations for each variant
    # annotated_vars = get_variant_information(vcf, vep_annotations, inclusion_variants, exclusion_variants, intervals)    
    
    ## Combine data
    my_list = []
    
    for var in annotated_vars: 
        my_list.append(var)

    ## Convert list of dictionaries to json file
    variant_path = args.var_path + '/gnomad_v3_variants.json'
    with open(variant_path, 'w') as fh:
        json.dump(annotated_vars, fh, indent=2)
    
#%%

if __name__ == '__main__':
    process_gnomad_v3_variants()
