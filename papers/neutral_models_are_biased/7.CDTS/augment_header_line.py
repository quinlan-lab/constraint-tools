import sys 

old_header_line = sys.stdin.readline().strip()
fields = old_header_line.split()
fields += ['feature', 'feature_chromosome', 'feature_start', 'feature_end', 'window_feature_overlap_bps']
new_header_line = '\t'.join(fields)
print(new_header_line)