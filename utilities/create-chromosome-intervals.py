# move to utils

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 12:59:04 2021

@author: u1240855
"""

import pandas as pd
import argparse

#%%
def get_chr_sizes(chr_sizes_file): 
    
    ## Initialize dictionary of chromosome sizes
    chr_sizes_dict = {}
    
    ## Save chr sizes as a dictionary
    with open(chr_sizes_file) as chr_sizes: 
        for line in chr_sizes: 
            (key, value) = line.split()
            chr_sizes_dict[str(key)] = int(value)
            
    return chr_sizes_dict
#%%     

#%%     
## Iterate through each chr size and get intervals
def split_chr_sizes(chr_sizes_dict, args):
    
    chr_sizes_split = []
    
    ## Initializze lists
    chromosome = []
    seg_start = []
    seg_end = []
    
    for chrom in chr_sizes_dict.keys():
        
        ## Get the total chromosome size
        chrom_size = int(chr_sizes_dict[chrom])

        ## Define segment length to interval the chrom sizes into user defined bins
        segment = int(chrom_size/int(args.bin_num))
        print("Bin number is: ", args.bin_num)
        
        ## Get the initial start/end to iteratively add to
        chromosome.append(chrom)
        seg_start.append(1)
        seg_end.append(seg_start[-1] + segment)
        
        
        while seg_end[-1] < chrom_size: ## While the interval's end is less than the total chr size
            chromosome.append(chrom)
            seg_start.append(seg_end[-1] + 1)
            seg_end.append(seg_end[-1] + segment + 1)
            
        ## Make the last end breakpoint the total chr length
        seg_end[-1] = chrom_size
                
    ## Convert chromosome intervals to a pandas dataframe
    chr_sizes_split = pd.DataFrame({'chromosome': chromosome, 'start': seg_start, 'end': seg_end})
    
    ## Output chromosome intervals
    output = args.output
    chr_sizes_split.to_csv(output, sep="\t", index=False)
    
#%% 

#%%
def parse_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--chr-sizes-file', type=str, dest='chr_sizes_file', help='')
    parser.add_argument('--bin-num', type=str, dest='bin_num', help='')
    parser.add_argument('--output', type=str, help='')
    
    return parser.parse_args()
#%%
 
if __name__ == '__main__': 
    
    ## Read in arguments
    args = parse_arguments()
    
    ## Store chr sizes in dictionary
    chr_sizes_dict = get_chr_sizes(args.chr_sizes_file)
    
    ## Interval the chr sizes
    split_chr_sizes(chr_sizes_dict, args)
