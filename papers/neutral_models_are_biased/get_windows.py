CONSTRAINT_TOOLS = '/scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools'
CONSTRAINT_TOOLS_DATA = '/scratch/ucgd/lustre-labs/quinlan/data-shared/constraint-tools'

import sys
sys.path.append(f'{CONSTRAINT_TOOLS}/utilities')

import pandas as pd
from functools import reduce
import numpy as np 

def compute_N_mean_null_gnocchi(row): 
    a = 1 
    b = -(2*row['N_observed'] + row['gnocchi']**2)
    c = row['N_observed']**2
    sqrt = np.sqrt(b**2 - 4*a*c)
    sign = 1 if row['gnocchi'] > 0 else -1
    return (-b + sign*sqrt)/(2*a)    

# Non-exonic windows, with Gnocchi and various features (e.g. GC content), and enhancer-overlap status 
def get_windows_with_GC_content_and_cpg_islands(): 
  print(CONSTRAINT_TOOLS_DATA)
  df1 = pd.read_csv(
    f'{CONSTRAINT_TOOLS_DATA}/chen-et-al-2023-published-version/41586_2023_6045_MOESM4_ESM/Supplementary_Data_2.gnocchi.N_expected.N_observed.B.paternal_recombination_rate.maternal_recombination_rate.gBGC-tract-counts.non-exonic.gBGC.bed', 
    sep='\t', 
  )
  df1 = df1.drop(columns=['N_expected']) 

  for gc_window_size in [
    1000, 
    10000, 
    100000,
    1000000,
  ]:
    df_temp = pd.read_csv(
      f'{CONSTRAINT_TOOLS_DATA}/chen-et-al-2023-published-version/41586_2023_6045_MOESM4_ESM/Supplementary_Data_2.GC_content_{gc_window_size}.bed', 
      sep='\t', 
    )
    df_temp = df_temp[['chen_chrom', 'chen_start', 'chen_end', 'window_GC_content']]
    df_temp = df_temp.rename(columns={
      'chen_chrom': 'chrom', 
      'chen_start': 'start', 
      'chen_end': 'end',
      'window_GC_content': f'GC_content_{gc_window_size}bp'
    })
    df1 = pd.merge(df1, df_temp, on=['chrom', 'start', 'end'], how='inner')

  # created using: experiments/germline-model/chen-et-al-2022/cpg-island-enrichment.ipynb
  df2 = pd.read_csv(
    f'{CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/mchale.kmerSizes.trainSets.noisy.enhancer-exon-cpgIsland.bed',
    sep='\t', 
  )
  df2 = df2[['chromosome', 'start', 'end', 'cpg_island overlap', 'window overlaps cpg_island']]
  df2 = df2.rename(columns={
    'chromosome': 'chrom', 
    'cpg_island overlap': 'cpg_island_overlap', 
    'window overlaps cpg_island': 'window_overlaps_cpg_island'
  })

  dfs = [df1, df2]
  df = reduce(lambda left, right: pd.merge(left, right, on=['chrom', 'start', 'end'], how='inner'), dfs)

  df['N_mean_null_gnocchi'] = df.apply(compute_N_mean_null_gnocchi, axis=1)

  return df

