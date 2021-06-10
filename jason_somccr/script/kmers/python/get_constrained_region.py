#!/usr/bin/env python3

'''
Python function to intervals of constraint from the variant file
'''

#%%
import pandas as pd

#%%

'''
Pseudocode: 
    1) Read in the variant file 
    2) Apply filters for genes/cancer types of interest
    3) Sort the variants by chr --> start --> stop
    4) Get the unique vales of the intervals
    5) Calculate intervals under constraint
'''

## Define the function
def get_constrained_region(nonsilent_df):
    
    ## Get columns of interest
    nonsilent_df = nonsilent_df[['chromosome', 'start', 'stop']]
    
    ## Apply filters here: 
    
        
    ## Sort dataset by chr --> start --> stop
    nonsilent_df_sorted = nonsilent_df.sort_values(by=['chromosome', 'start', 'stop'])
    
    ## Get unique variant positions
    nonsilent_df_sorted_uniq = nonsilent_df_sorted.drop_duplicates()
    
    ## Initialize list to store variant information
    var_list = []
    
    ## Iterate through each unique variant position and get the 
    for chromosome in range(1,25): ## Iterate through each chromosome
    
        ## Print status message
        print(f'Getting constrained intervals found on chromosome{chromosome}')
        
        ## Subset chromosome data
        chr_df = nonsilent_df_sorted_uniq.loc[(nonsilent_df_sorted_uniq['chromosome'] == chromosome)]
        
        '''
        Intervals NEED to be 1-based coordinates to get the nucleotide sequence
        
        Variant 1 = 1:321-322
        Variant 2 = 1:348-349
        
        Stop of v1 is in 1-based coordinates --> 322 (1-based) = 322 (1-based)
        Start of v2 is in 0-based coordinates --> 348 (0-based) = 348 (1-based)
        '''
        #chr_df = chr_df.head()
        
        ## Get all the interval start positions
        interval_start = chr_df['stop'].tolist()
        interval_start.pop()
        
        ## Get all the inrerval stop positions
        interval_stop = chr_df['start'].tolist()
        interval_stop.pop(0)
        
        ## Get number of variants
        num_var = chr_df.shape[0] -1
        
        ## Define intervals to calculate kmers from the constrained interval
        for num in range(0, num_var):
            
            ## Define the interval's region for nucleotide sequence
            chrom = int(chromosome)
            start = interval_start[num] - 1
            stop = interval_stop[num] + 1
            
            ## Make list
            var_info = [chrom, start, stop]
            var_list.append(var_info)
        
    ## Store variant information into a pandas dataframe
    chr_var_df = pd.DataFrame(var_list, columns = ['chromosome', 'start', 'stop'])
    
    ## Print status message
    print('Writing the output file')
    
    ## Write the data to an output file
    output_file = '/scratch/ucgd/lustre-work/quinlan/u1240855/somccr/data/output/constrained_intervals/constrained_interval.bed'
    chr_var_df.to_csv(output_file, sep="\t", header=True, index=False)   
    
#%%

## Read in the cds bed coordinates as a pandas dataframe
variant_file = '/scratch/ucgd/lustre-work/quinlan/u1240855/somccr/data/output/tcga_icgc.sorted.filtered.snv.CDS.knownCDSlength.no_chr_prefix.fix_header.renameXY.bed'
variant_df = pd.read_csv(variant_file, sep="\t", usecols=('chromosome', 'start', 'stop', 'variant_class', 'alt'))

## Select for non-silent mutations
nonsilent_df = variant_df.loc[(variant_df['variant_class'] != 'Silent') & (variant_df['alt'].isin(['A', 'C', 'G', 'T']))]

## Get the intervals of constraint from the variant file 
get_constrained_region(nonsilent_df = nonsilent_df)