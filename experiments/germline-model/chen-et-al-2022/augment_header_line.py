old_header_line = 'chromosome start end     N_observed      K_observed      M       N_bar_3 K_bar_3 N_bar_5 K_bar_5 N_bar_7 K_bar_7 chen_zscore'
fields = old_header_line.split()
fields += ['feature', 'feature_chromosome', 'feature_start', 'feature_end', 'window_feature_overlap_bps']
new_header_line = '\t'.join(fields)
print(new_header_line)