CONSTRAINT_TOOLS = '/scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools'
CONSTRAINT_TOOLS_DATA = '/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools'

import pandas as pd 
import warnings

pd.set_option('display.max_columns', 50)

NUMBER_WINDOWS_PER_QUANTILE = 10

def get_windows_scores_annotations():
    filename = f'{CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/enhancer-characteristics-enrichment.bed'
    df = pd.read_csv(filename, sep='\t')
    df = df[df['window overlaps merged_exon'] == False]
    return df

windows_scores_annotations_noncoding = get_windows_scores_annotations()

def compute_array_of_quantiles():
  array_of_quantiles = [0.00, 0.005, 0.01, 0.02, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.98, 0.99, 0.995, 1.00]
  # array_of_quantiles = [0.00, 0.01, 0.02, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.98, 0.99, 1.00]
  starts = array_of_quantiles[:-1]
  ends = array_of_quantiles[1:]
  quantile_labels = [f'{start} - {end}' for start, end in zip(starts, ends)]
  return array_of_quantiles, quantile_labels

def label_windows_with_score_quantiles_core(df, score): 
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")        
    array_of_quantiles, quantile_labels = compute_array_of_quantiles()
    df[f'{score} quantile'], bins = pd.qcut(
      df[score],
      q = array_of_quantiles, 
      labels = quantile_labels, 
      retbins = True,
#         duplicates='drop'
    )
    return df

def subset(df, score='negative new chen zscore'):
  df = label_windows_with_score_quantiles_core(df, score)
  df = (
    df
    .groupby([f'{score} quantile'], as_index=True)
    .apply(lambda df_subset: df_subset.sample(NUMBER_WINDOWS_PER_QUANTILE))
    .reset_index(drop=True)
  )
  out_filename = f'{CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/enhancer-characteristics-enrichment-subset.bed'
  df.to_csv(                                                                                
      out_filename,
      sep = '\t',
      index = False
  )
  print(f'Wrote data to {out_filename}')

subset(windows_scores_annotations_noncoding)

