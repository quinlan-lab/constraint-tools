import sys 

input_filename = sys.argv[1] 
realization = sys.argv[2] 

import pandas as pd 

chen_mchale_DELs = pd.read_csv(input_filename, sep='\t')

chen_mchale_DELs['negative chen zscore'] = -chen_mchale_DELs['chen_zscore']

def compute_window_overlaps_DEL(df):     
  for source in ['CCDG', 'gnomAD', '1000G']:
    df[f'window overlaps {source} DEL'] = df[f'{source}_DELs'] > 0
  return df 

chen_mchale_DELs = compute_window_overlaps_DEL(chen_mchale_DELs)

def rename(df):
  for source in ['CCDG', 'gnomAD', '1000G']:
    old = f'{source}_bp_overlap'
    new = f'fraction of window that intersects {source} DEL'
    df = df.rename(columns={old: new})
  return df

chen_mchale_DELs = rename(chen_mchale_DELs)

import warnings

def label_windows_with_score_quantiles(df, score): 
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")        

    array_of_quantiles = [0.00, 0.01, 0.02, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.98, 0.99, 1.00]
    starts = array_of_quantiles[:-1]
    ends = array_of_quantiles[1:]
    quantile_labels = [f'{start} - {end}' for start, end in zip(starts, ends)]

    df[f'{score} quantile'], bins = pd.qcut(
      df[score],
      q = array_of_quantiles, 
      labels = quantile_labels, 
      retbins = True,
#         duplicates='drop'
    )

    return df

chen_mchale_DELs = label_windows_with_score_quantiles(
  chen_mchale_DELs, 
  score='N_bar'
)
chen_mchale_DELs = label_windows_with_score_quantiles(
  chen_mchale_DELs, 
  score='negative chen zscore'
)

def aggregate(df, zscore, source):
  groups = df.groupby([f'{zscore} quantile'])
  aggregated = groups.agg({
    f'window overlaps {source} DEL': ['mean'],
  })
  return aggregated

pd.set_option('display.max_columns', 30)

import numpy as np 

def save_window_fractions_grouped_by(zscore):
  dfs = []
  for source in ['CCDG', 'gnomAD', '1000G']:
    dfs.append(aggregate(
      df = chen_mchale_DELs, 
      zscore = zscore,
      source = source
    ))
    
  from functools import reduce
  df = reduce(lambda x, y: pd.merge(x, y, on=f'{zscore} quantile'), dfs)

  output_filename = f'fraction_windows_that_overlap_DEL__grouped_by_{zscore}_quantile__realization_{realization}.pkl'
  df.to_pickle(output_filename)
  # print(pd.read_pickle(output_filename))

save_window_fractions_grouped_by('N_bar')
save_window_fractions_grouped_by('negative chen zscore')

