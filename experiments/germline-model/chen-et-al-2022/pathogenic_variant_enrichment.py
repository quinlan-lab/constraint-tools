import pandas as pd 
import numpy as np 
import warnings
from collections import defaultdict
from tqdm import tqdm 
import sys 
import pickle

pd.set_option('display.max_columns', 30)

CONSTRAINT_TOOLS_DATA = '/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools'
THRESHOLD_CLINVAR_COUNT = 4 

def get_chen_mchale_clinvar():
  filename = f'{CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/chen-mchale-enhancer-exon-clinvar.bed'
  chen_mchale_clinvar = pd.read_csv(filename, sep='\t')
  return chen_mchale_clinvar

def get_chen_mchale_mendelian_variants():
  filename = f'{CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/chen-mchale-enhancer-exon-mendelian-variants.bed'
  chen_mchale_mendelian_variants = pd.read_csv(filename, sep='\t')
  return chen_mchale_mendelian_variants

def wrangle_K_bar(df): 
    df['K_bar'] = pd.to_numeric(df['K_bar'], errors='coerce') # convert '.' to 'NaN'
    #     df = df.dropna(subset=['K_bar']) # drop windows for which K_bar == NaN
    df['negative K_bar'] = -df['K_bar']
    return df
    
def get_chen_mchale_clinvar_mendelianVariants():
  chen_mchale_clinvar = get_chen_mchale_clinvar()
  chen_mchale_mendelian_variants = get_chen_mchale_mendelian_variants()
  chen_mchale_clinvar_mendelianVariants = chen_mchale_clinvar.merge(chen_mchale_mendelian_variants)
  chen_mchale_clinvar_mendelianVariants = wrangle_K_bar(chen_mchale_clinvar_mendelianVariants)
  return chen_mchale_clinvar_mendelianVariants

def label_windows_with_score_quantiles_core(df, score): 
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

def label_windows_with_score_quantiles(df):   
  df = label_windows_with_score_quantiles_core(df, score='N_bar')
  df = label_windows_with_score_quantiles_core(df, score='negative K_bar')
  df = label_windows_with_score_quantiles_core(df, score='negative chen zscore')
  return df 

def aggregate(df, zscore):
    groups = df.groupby([f'{zscore} quantile'])
    aggregated = groups.agg({
        'ClinVar SNV count': ['mean', 'std', 'count', 'sum'],
        'Mendelian variant count': ['mean', 'std', 'count', 'sum'],
        'merged_exon overlap': ['mean', 'std', 'count']
    })
    return aggregated

def resample(df): 
    return df.sample(n=len(df), replace=True)

def compute_enrichment(df, zscore, variant_type, do_resample): 
  if do_resample: df = resample(df)
  aggregated = aggregate(df, zscore)
  zscore_quantiles = np.array(aggregated.index)    
  window_counts = np.array(aggregated[(f'{variant_type} count', 'count')])
  variant_counts_observed = np.array(aggregated[(f'{variant_type} count', 'sum')])
  variant_counts_total = np.sum(variant_counts_observed)
  zscore_quantile_probabilities_null = window_counts/np.sum(window_counts) 
  mean_variant_counts_null = variant_counts_total*zscore_quantile_probabilities_null
  enrichment = variant_counts_observed/mean_variant_counts_null
  return zscore_quantiles, enrichment 

def bootstrap_enrichment(df, zscore, variant_type, number_bootstrap_samples): 
  samples = []
  for i in tqdm(range(number_bootstrap_samples), desc=f'{zscore}; {variant_type}'):
    zscore_quantiles, enrichment = compute_enrichment(df, zscore, variant_type, do_resample=True)
    samples.append(enrichment)
  samples = np.vstack(samples)
  ci_lower = np.percentile(samples, 2.5, axis=0) 
  mean = np.mean(samples, axis=0) 
  ci_upper = np.percentile(samples, 97.5, axis=0) 
  return {
    'zscore_quantiles': zscore_quantiles,
    'mean': mean,
    'ci_lower': ci_lower, 
    'ci_upper': ci_upper,    
  }
 
def show_outliers(df): 
    df = df.copy()
    df = df[
        (df['ClinVar SNV count'] > THRESHOLD_CLINVAR_COUNT) 
    ]
    df['region'] = df['chromosome'] + ':' + df['start'].apply(str) + '-' + df['end'].apply(str)
    return df

def remove_outliers(df):
    df = df.copy()
    df = df[df['ClinVar SNV count'] < THRESHOLD_CLINVAR_COUNT]
    return df 

def save_bootstrapped_enrichment(number_bootstrap_samples): 
  chen_mchale_clinvar_mendelianVariants = get_chen_mchale_clinvar_mendelianVariants()
  chen_mchale_clinvar_mendelianVariants = label_windows_with_score_quantiles(chen_mchale_clinvar_mendelianVariants)
  chen_mchale_clinvar_mendelianVariants = remove_outliers(chen_mchale_clinvar_mendelianVariants)

  bootstrapped_enrichment = defaultdict(lambda: {})
  for zscore in ['N_bar', 'negative K_bar', 'negative chen zscore']: 
    for variant_type in ['ClinVar SNV', 'Mendelian variant']:
      bootstrapped_enrichment[variant_type][zscore] = bootstrap_enrichment(
        chen_mchale_clinvar_mendelianVariants, 
        zscore, 
        variant_type, 
        number_bootstrap_samples
      )
  bootstrapped_enrichment = dict(bootstrapped_enrichment)
  with open('pathogenic_variant_enrichment.pkl', 'wb') as handle:
      pickle.dump(bootstrapped_enrichment, handle, protocol=pickle.HIGHEST_PROTOCOL)

if __name__ == '__main__':
  save_bootstrapped_enrichment(number_bootstrap_samples=int(sys.argv[1]))

