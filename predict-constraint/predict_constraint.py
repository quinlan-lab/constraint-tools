#!/usr/bin/env python3

'''
Steps to predict constraint for a given gene
1) Define gene of interest
2) Query gtf file for gene of interest
3) identify coordinates for coding exons
4) obtain sequence inforation for each exon --> flatten into a meta exon
5) Tile along meta exon 
    - compute observed/expected
    - get difference between observed/expected
'''

#%%
import pysam
import pandas as pd
import numpy as np
import json
import argparse
import re
from scipy.stats import poisson

from pack_unpack import pack, unpack
from compute_mutation_counts import create_windows, compute_expected_mutation_count
from kmer import fetch_kmer_from_sequence, middle_base, middle_index, check_for_Ns
from fetch_SNVs import fetch_SNVs
from colorize import print_json, print_string_as_info, print_string_as_info_dim, print_unbuffered, print_string_as_error


#%%

#%%
def parse_arguments(): 
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--transcript', type=str, dest='transcript', help='')
    parser.add_argument('--transcript-intervals', type=str, dest='transcript_intervals', help='')
    parser.add_argument('--window-size', type=int, default=101, dest='window_size', help='')
    parser.add_argument('--window-stride', type=int, default=50, dest='window_stride', help='')
    parser.add_argument('--model-filename', type=str, dest='model_filename', help='')
    parser.add_argument('--output-path', type=str, dest='output_path', help='')
    
    return parser.parse_args()
#%%

#%%
# Test input

## Define gene of interest
transcript = "ENST00000011619"

## Define transcript intervals
transcript_intervals="/scratch/ucgd/lustre-work/quinlan/u1240855/constraint-tools/predict-constraint/intermediate_files/intervals/ENST00000011619/ENST00000011619.sorted.final.bed"
intervals = pd.read_csv(transcript_intervals, sep="\t", header=None)
intervals.columns = ['chrom', 'start', 'end', 'tag']

## Read in the json model file
model_filename = "/scratch/ucgd/lustre-work/quinlan/u1240855/constraint-tools/tests/model.json"

## Define output directory
output_path = '/scratch/ucgd/lustre-work/quinlan/u1240855/constraint-tools/predict-constraint/intermediate_files/transcript_predict_constraint'

## Define window parameters
window_size = 101
window_stride = 50
#%%

#%%
def get_cDNA_coord(pos, seq): 
    ## Get the length of the interval sequence
    seq_length = len(seq)
    
    start_pos = pos
    end_pos = pos + seq_length
    
    return start_pos, end_pos
#%%

#%%
def get_interval_sequences(intervals, genome, transcript): 
    
    ## Initialize dictionary
    interval_dict = {}
    regions = {}
    
    ## Initialize length of transcript
    transcript_length = 0
    transcript_length_covered = 0
    transcript_length_notcovered = 0
    
    ## Initialize the cDNA start position
    cDNA_pos = 0
    
    ## Initialize final sequence
    final_seq = ""
    
    ## Iterate through each interval and get the sequence
    for index, interval in intervals.iterrows(): 
        
        ## Define the region
        region = pack(interval['chrom'], interval['start'], interval['end'])
        print(region)
        
        ## Get the sequence of the interval
        seq = check_for_Ns(genome.fetch(interval['chrom'], interval['start'], interval['end']))

        ## Check to see if the seq is covered 
        if interval['tag'] != " covered": 
            
            ## Remove nucleotide information
            print("Interval {} is not covered... Removing the sequence information...".format(region))
            seq = '-' * len(seq)
            
            ## Add to not covered transcript length
            transcript_length_notcovered += len(seq)
            
        if interval['tag'] == " covered": 
            ## Add to covered transcript length
            transcript_length_covered += len(seq)
            
        ## Get the cDNA position for the interval
        start_window, end_window = get_cDNA_coord(pos = cDNA_pos, seq = seq)

        ## Store sequence and cDNA window information 
        sub_dict = {'seq': seq,
                    'start': start_window,
                    'end': end_window}
        
        ## Update the cDNA position
        cDNA_pos = end_window 
        
        ## Store all information within each region
        regions[region] = sub_dict
        
        ## Add to transcript length
        transcript_length += len(seq)
        
        ## Combine sequences
        final_seq = final_seq + seq

    ## Fill in the dictionary
    interval_dict['transcript_length'] = transcript_length
    interval_dict['transcript_length_covered'] = transcript_length_covered
    interval_dict['transcript_length_notcovered'] = transcript_length_notcovered
    
    ## Combine the sequences from all intervals
    interval_dict['chromosome'] = interval['chrom']
    interval_dict['transcript'] = transcript
    interval_dict['transcript_seq'] = final_seq
    
    ## Add sequence information for each interval
    interval_dict['regions'] = regions
    
    print("The transcript {0} has a length of {1} where {2} sites are covered and {3} sites are not covered... Making predictions on the covered regions...".format(transcript, transcript_length, transcript_length_covered, transcript_length_notcovered))
    
    return interval_dict
#%%

#%%
def map_variants_to_cDNA(interval_dict, mutations, genome, model): 
    
    ## Initialize list of mutations
    mutations_list = []
    
    ## Iterate through each region in the interval dictionary --> get mutations within that genomic space
    for region, value in interval_dict['regions'].items(): 
        region_mutations = fetch_SNVs(mutations, genome, region, meta=model, return_coords=True)
        
        ## Get the start/end gDNA coordinates for the region
        chromosome, region_start, region_end = unpack(region)
        
        ## Get the start/end cDNA coordinates for the region
        region_cDNA_start = interval_dict['regions'][region]['start']
        region_cDNA_end = interval_dict['regions'][region]['end']
        
        ## Map each mutation to cDNA coordinates
        if len(region_mutations) == 0: 
            print('No mutations found within the region: {}'.format(region))
        
        for mut in region_mutations: 
            ## Get the genomic start coordinate
            gDNA_start = int(mut['start'])
            
            ## Get cDNA offset
            cDNA_offset = gDNA_start - region_start
            
            ## Perform the map
            cDNA_start = region_cDNA_start + cDNA_offset
            cDNA_end = cDNA_start + 1
            
            ## Add cDNA coordinates to mutation dictionary
            mut['cDNA_start'] = cDNA_start
            mut['cDNA_end'] = cDNA_end
            mut['gDNA_region'] = region
            mut['cDNA_region'] = chromosome + ':' + str(region_cDNA_start) + "-" + str(region_cDNA_end)
            
            ## Append to mutations list
            mutations_list.append(mut)
            
    ## Convert to pandas dataframe
    final_mut = pd.DataFrame(mutations_list)
    
    return final_mut
#%%

#%%
def create_cDNA_windows(interval_dict, window_size, window_stride): 
    
    ## Define the region as the cDNA space
    region = interval_dict['chromosome'] + ":" + '0' + "-" + str(interval_dict['transcript_length'])
    
    ## Generate windows in cDNA space
    cDNA_windows = create_windows(window_size, window_stride, region)
    
    return cDNA_windows
#%%

#%%
def get_expected_mutation_count(cDNA_windows, interval_dict, model): 

    ## Initialize windows
    windows = []
    
    ## Iterate through each window and calculate the number of expected mutations
    for window in cDNA_windows: 
        
        ## Get the cDNA window start and stop
        window_chr = window['region']['Chromosome']
        window_start = window['region']['Start']
        window_end = window['region']['End']
        print("Computing expected mutation count for the cDNA window {}:{}-{}".format(window_chr, window_start, window_end))
        
        ## Get the entire sequence of the window
        seq = interval_dict['transcript_seq']
        
        ## Initialize expected mutation count for the window
        window_expected_mutation_count = 0
        
        ## Iterate through each position in the region
        for position in np.arange(window_start, window_end, 1):
            
            print(position)
            try:
                ## Get the position's kmer
                kmer = fetch_kmer_from_sequence(seq, position, model['kmer_size'])
                
                if len(kmer) < model['kmer_size']: 
                    print('Removing kmers with length less than the designated kmer size')
                    continue
                
            ## Check to see if a kmer could be determined
            except IndexError: 
                print_string_as_info_dim('IndexError at position: {}... \
                                         Could not determine kmer due to lack of flanking cDNA nucleotides...'.format(position))
            #if kmer == '': 
            #    print("Could not determine kmer due to lack of flanking cDNA nucleotides...")
            #    continue
                
            ## Check to see if the entire kmer sequence is covered
            if re.search('-', kmer): 
                print("At least 1+ site of the kmer {} is insufficiently covered...".format(kmer))
                continue
            
            ## Compute the mutation count for the kmer
            window_expected_mutation_count += compute_expected_mutation_count(kmer, model)

        
        ## Store expected mutation count for the given interval
        window['expected_mutation_count'] = window_expected_mutation_count
        
        ## Append window to windows list
        windows.append(window)
        
    return windows
#%%

#%%
def get_observed_mutation_count(cDNA_windows_exp, model, cDNA_variants): 
    ## Initialize windows
    windows = []
    
    ## Itearte through each window and calculate the number of observed mutations
    for window in cDNA_windows_exp:
        
        ## Get the window information
        window_chr = window['region']['Chromosome']
        window_start = window['region']['Start']
        window_end = window['region']['End']
        
        ## Filter cDNA variants to determine number of mutations within the window
        window_observed_mutation_count = cDNA_variants[(cDNA_variants['cDNA_start'] > window_start) & 
                      (cDNA_variants['cDNA_end'] < window_end)].shape[0]

        ## Store the observed mutation count value
        window['observed_mutation_count'] = window_observed_mutation_count
        
        ## Append window data to windows list
        windows.append(window)
    
    ## Return the data
    return windows
#%%

#%%
def compute_observed_expected_ratio(cDNA_windows_exp_obs): 
    
    ## Initialize windows
    windows = []
    
    ## Iterate through each window and calculate observed-expected ratio
    for window in cDNA_windows_exp_obs: 
                
        ## Compute the window size
        chromosome, start, end = window['region']['Chromosome'], window['region']['Start'], window['region']['End']
        window['chrom'] = chromosome
        window['start'] = start
        window['end'] = end
        window['length'] = end - start
    
        ## Calculate poisson probability
        prob = poisson.pmf(k = window['observed_mutation_count'], mu = window['expected_mutation_count'])
        window['probability'] = prob
        
        ## Check to see if expected value is 0
        if window['expected_mutation_count'] == 0: 
            window['o_e_ratio'] = 0
            windows.append(window)
            continue
        
        ## Determine the observed-expected ration
        window['o_e_ratio'] = window['observed_mutation_count']/window['expected_mutation_count']
                
        ## Append data
        windows.append(window)
    
    ## return data
    return windows
#%%

#%%
def cDNA_to_genomic(interval_dict, transcript_windows_final): 
     ## Reorganize the interval dictionary data
    interval_df = []
    
    for interval, value in interval_dict['regions'].items(): 
        
        ## Initialize coverage tag
        cov_status = "covered"
        
        if re.search("-", interval_dict['regions'][interval]['seq']): 
            cov_status = "not_covered"
        
        ## Get the genomic coordinates of the interval
        g_chr, g_start, g_end = unpack(interval)
        
        ## Get the cDNA start/end coordinates
        c_start = value['start']
        c_end = value['end']
        
        ## Append data        
        interval_df.append({'Chromosome': g_chr, 'Genomic_start': g_start, 'Genomic_end': g_end,
         'cDNA_start': c_start, 'cDNA_end': c_end, 'cov_status': cov_status})
        
        
    return pd.DataFrame(interval_df)
#%%

#%%
def cDNA_to_genomic_offset(subset_mapped_intervals, window_start, window_end): 
    ## Get the overhanging offset for the start position
    left_pos = subset_mapped_intervals[(subset_mapped_intervals['cDNA_start'] <= window_start) & 
                                       (subset_mapped_intervals['cDNA_end'] >= window_start)]
    left_offset = window_start - left_pos['cDNA_start']
    new_genomic_start = left_pos['Genomic_start'] + left_offset
    
    ## Get the overhanging offset for the end position
    right_pos = subset_mapped_intervals[(subset_mapped_intervals['cDNA_start'] <= window_end) &
                                        (subset_mapped_intervals['cDNA_end'] >= window_end)]
    right_offset = right_pos['cDNA_end'] - window_end
    new_genomic_end = right_pos['Genomic_end'] - right_offset
    
    return new_genomic_start, new_genomic_end
#%%

#%%
def map_windows_to_genomic(interval_dict, transcript_windows_final, transcript): 
        
    ## Map cDNA intervals to genomic coordinates
    mapped_intervals = cDNA_to_genomic(interval_dict, transcript_windows_final)
    
    ## Initialize windows
    df = pd.DataFrame()
    
    ## Iterate through each window --> get overlapping genomic coordinates
    for window in transcript_windows_final: 
        print(window)
        print(window['probability'])
        
        ## Get the cDNA coordinates for the window
        window_chr, window_start, window_end = window['chrom'], window['start'], window['end']
        
        ## Remove intervals not sufficiently covered
        # covered_intervals = mapped_intervals[mapped_intervals['cov_status'] == "covered"]
        
        ## Get the overlapping genomic intervals with the cDNA window
        subset_mapped_intervals = mapped_intervals[((mapped_intervals['cDNA_start'] <= window_start) & (mapped_intervals['cDNA_end'] >= window_start )) | 
                                                   ((mapped_intervals['cDNA_start'] >= window_start) & (mapped_intervals['cDNA_end'] <= window_end)) | 
                                                   ((mapped_intervals['cDNA_start'] <= window_end) & (mapped_intervals['cDNA_end'] >= window_end))]
        
        ## Change genomic start/end if there is overhanging genomic space relative to the cDNA window
        new_genomic_start, new_genomic_end = cDNA_to_genomic_offset(subset_mapped_intervals, window_start, window_end)
        
        ## Rename index names
        row_num = subset_mapped_intervals.shape[0]
        subset_mapped_intervals.index = range(row_num)
        
        ## Update the mapped intervals 
        subset_mapped_intervals.at[0, 'Genomic_start'] = new_genomic_start
        subset_mapped_intervals.at[0, 'cDNA_start'] = window_start
        subset_mapped_intervals.at[row_num - 1, 'Genomic_end'] = new_genomic_end
        subset_mapped_intervals.at[row_num - 1, 'cDNA_end'] = window_end
        
        ## Add window id and values to mapped genomic space
        subset_mapped_intervals['transcript'] = transcript
        subset_mapped_intervals['cDNA_region'] = pack(window['chrom'], window['start'], window['end'])
        subset_mapped_intervals['cDNA_window_obs'] = window['observed_mutation_count']
        subset_mapped_intervals['cDNA_window_exp'] = window['expected_mutation_count']
        subset_mapped_intervals['cDNA_window_obs_exp_ratio'] = window['o_e_ratio']
        subset_mapped_intervals['cDNA_window_length'] = window['length']
        subset_mapped_intervals['probability'] = window['probability']
        
        ## Append to window's pandas dataframe
        df = df.append(subset_mapped_intervals, ignore_index=True)
        
    return df

#%%

#%%
def write_output(transcript_windows_final, output_path, transcript): 
    fh = output_path + '/' + transcript + '_predict_constraint.bed' ## Verify bed coordinates
    
    ## Write to output path
    df.to_csv(fh, sep='\t', header=True, index=False)
#%%

#%%
## Define function to predict constraint
def predict_constraint(): 

    ##################################    
    ##### Define input variables #####
    ##################################
    
    ## Get the arguments
    args = parse_arguments()
    
    ## Define transcript of interest
    transcript = args.transcript
    print("Predicting constrained regions for transcript: {}".format(transcript))
    
    ## Define window size and stride
    window_size = args.window_size
    window_stride = args.window_stride
    
    ## Read in the intervals data
    intervals = pd.read_csv(args.transcript_intervals, sep="\t", header=None)
    intervals.columns = ['chrom', 'start', 'end', 'tag']
    
    ## Read in data from the model 
    with open(args.model_filename) as fh: 
        model = json.load(fh)
    mutations = pysam.TabixFile(model['mutations'])
    genome = pysam.FastaFile(model['genome'])
    
    ##############################:
    ##### Predict constraint #####
    ##############################
    
    ## Get the sequence of each covered/notcovered interval for the transcript
    print("Identifying sequence of covered and uncovered sites for the transcript {}...".format(transcript))
    interval_dict = get_interval_sequences(intervals, genome, transcript)
    
    ## Map variants to cDNA coordinates
    print("Mapping variants lying within coding exons of {} to cDNA coordinates...".format(transcript))
    cDNA_variants = map_variants_to_cDNA(interval_dict, mutations, genome, model)
    
    ## Create windows to predict constraint on using cDNA coordinates
    print("Generating windows using cDNA coordinates to predict constraint on...")
    cDNA_windows = create_cDNA_windows(interval_dict, window_size, window_stride)
    
    ## Get the expected mutation count for each interval in the transcript
    print("Computing the expected number of mutations within each cDNA interval...")
    cDNA_windows_exp = get_expected_mutation_count(cDNA_windows, interval_dict, model)
    
    ## Get the observed mutation count for each interval in the transcript
    print("Computing the observed number of mutations within each cDNA interval...")
    cDNA_windows_exp_obs = get_observed_mutation_count(cDNA_windows_exp, model, cDNA_variants)
    
    ## Create windows to calculate o/e using cDNA coordinates
    print("Computing the o/e ratio within each cDNA interval...")
    transcript_windows_final = compute_observed_expected_ratio(cDNA_windows_exp_obs)
    
    ## Map cDNA windows to genomic coordinates
    print("Mapping cDNA windows to genomic coordinates...")
    df = map_windows_to_genomic(interval_dict, transcript_windows_final, transcript)
    
    ## Write the transcript_windows data to output_path
    write_output(df, output_path, transcript)
    
    ## Print success signature
    print("Predicted constrained regions for transcript: {}".format(transcript))
    print("Success!")

#%%

if __name__ == '__main__': 
    predict_constraint()
    
    
#%%
for 
#%%


