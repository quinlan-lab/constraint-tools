#!/usr/bin/env python3

"""
Function to download chromsome fasta files (hg19)
"""

from os import path
import concurrent.futures
import subprocess

def get_chr_fasta_file(chromosome):
        
    ## Define the file name
    filename = "/scratch/ucgd/lustre-work/quinlan/u1240855/somccr/data/reference/chr/chr" + chromosome + ".fa" 
    filename_index = "/scratch/ucgd/lustre-work/quinlan/u1240855/somccr/data/reference/chr/chr" + chromosome + ".fa.fai" 
    
    ## Check if the chr fasta file exists
    if path.exists(filename) == True and path.exists(filename_index) == True: 
        return f'{filename} already exists, moving on...'
        
    else:         
        ## Define the URL
        url = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr' + chromosome + '.fa.gz'
        
        ## Run the bash script to download UCSC fasta file (and gunzip)
        subprocess.run(['~/git/somccr/script/kmers/slurm/get_chr_fasta_files.sh %s %s' %(url, filename)], shell=True)
        
        return f'Downloaded UCSC fasta file for chr{chromosome}...'
    
with concurrent.futures.ThreadPoolExecutor() as executor:
    
    ## Define list of chromosomes
    my_list = list(range(1,22+1))
    my_list.extend(['X', 'Y'])
    my_list = list(map(str, my_list))
    
    ## Go through each chromosome and get the fasta file if needed
    executor.map(get_chr_fasta_file, my_list)
    




