import pandas as pd 

def get_noncoding_windows(): 
    filename = '/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/benchmark-genome-wide-predictions/chen-et-al-2022/chen-mchale-enhancer-exon.bed'
    df = pd.read_csv(filename, sep='\t')
    df = df[df['window overlaps merged_exon'] == False] 
    columns_to_drop = [
        'enhancer overlap',
        'merged_exon overlap', 
        'window overlaps enhancer',
        'window overlaps merged_exon',
        'window overlaps (enhancer, merged_exon)', 
        'negative chen zscore',
        'singleton proportion', 
        'singleton proportion error (ind. sites)',
        'negative N_bar'
    ]
    df = df.drop(columns = columns_to_drop)
    df['chromosome'] = df['chromosome'].apply(lambda chromosome: chromosome.strip('chr'))
    return df

import sys 
import errno

try: 
    get_noncoding_windows().to_csv(sys.stdout, sep='\t', header=False, index=False)
except IOError as e:
    if e.errno == errno.EPIPE: pass
    else: raise

