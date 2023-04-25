old_header_line = 'chromosome      start   end     N_observed      N_bar_7_noncoding       N_bar_1_noncoding       N_bar_1_noncoding_noisy N_bar_3_noncoding       N_bar_3_noncoding_noisy N_bar_5_noncoding       N_bar_5_noncoding_noisy'
fields = old_header_line.split()
fields += ['feature', 'feature_chromosome', 'feature_start', 'feature_end', 'window_feature_overlap_bps']
new_header_line = '\t'.join(fields)
print(new_header_line)