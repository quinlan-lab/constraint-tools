import pandas as pd 
import numpy as np 
from scipy.signal import savgol_filter
import warnings 
from collections import defaultdict
import pickle
import sys 
from tqdm import tqdm 

parameters_to_estimate_singleton_probability_under_null = { 
    'gnomAD': dict(
      number_length_quantiles=200, 
      savgol_window_length=51, 
      degree=3,
      array_of_quantiles=[0.000, 0.02, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.98, 1.00]
    ),
    'CCDG': dict(
      number_length_quantiles=100, 
      savgol_window_length=11, 
      degree=5,
      array_of_quantiles=[0.000, 0.02, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.98, 1.00]
    ),
    '1000G': dict(
      number_length_quantiles=100, 
      savgol_window_length=51, 
      degree=3,
      array_of_quantiles=[0.000, 0.02, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.98, 1.00]
    )
}

def get_noncoding_svs_windows(source):
    CONSTRAINT_TOOLS_DATA = '/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools'
    # created using ${CONSTRAINT_TOOLS}/experiments/germline-model/chen-et-al-2022/intersect-noncoding-svs-with-windows.sh : 
    filename = f'{CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/{source}-noncoding-svs-chen-mchale.kmerSizes.trainSets.enhancer-exon.bed'
    df = pd.read_csv(filename, sep='\t')
    df['negative new chen zscore'] = -df['new chen zscore']
    
    columns = ['sv_id', 'sv_length', 'alt_allele_count', 'negative new chen zscore']
    for kmer_size in [3, 5, 7]: 
        for train_set_label in ['coding', 'noncoding', 'chenWindows']: 
            columns.append(f'N_bar_{kmer_size}_{train_set_label}')            
    df = df[columns]
    
    df['sv is singleton'] = (df['alt_allele_count'] == 1) | (df['alt_allele_count'] == '1')  
    return df 

# zscore_percentile = 10 

# def mean_of_bottom(xs): 
#     xs = np.array(xs)
#     threshold = np.percentile(xs, zscore_percentile) 
#     bottom = xs[xs <= threshold] 
#     return np.mean(bottom)

def aggregate_over_windows(df): 
    groups = df.groupby(['sv_id', 'sv_length', 'alt_allele_count', 'sv is singleton'])
#     custom_aggregation_function = (f'mean of bottom {zscore_percentile}%', mean_of_bottom)

    aggregation_functions = {'negative new chen zscore': ['min', 'mean']}
    for kmer_size in [3, 5, 7]: 
        for train_set_label in ['coding', 'noncoding', 'chenWindows']: 
            aggregation_functions[f'N_bar_{kmer_size}_{train_set_label}'] = ['min', 'mean']
    aggregated = groups.agg(aggregation_functions)
    
    df = aggregated.reset_index()
    df.columns = [' '.join(col[::-1]).strip() for col in df.columns.values]
    return df

def resample(df): 
    return df.sample(n=len(df), replace=True)

def get_svs(source, do_resample): 
    df = get_noncoding_svs_windows(source)
    df = aggregate_over_windows(df)
    if do_resample: df = resample(df)
    return df
   
def label_svs_with_length_quantiles(df, params): 
    df['sv length quantile'] = pd.qcut(
        df['sv_length'],
        q = params['number_length_quantiles'], 
        labels = False, 
        retbins = False,
#         duplicates='drop'
    )
    return df

def aggregate_over_length_quantiles(df, params): 
    groups = df.groupby(['sv length quantile'])
    aggregated = groups.agg({
        'sv_length': ['mean', 'count'],
        'sv is singleton': ['mean', 'count']
    })
    aggregated = aggregated.reset_index() 
    aggregated['estimated singleton probability under null'] = savgol_filter(
        aggregated[('sv is singleton', 'mean')], 
        window_length = params['savgol_window_length'], 
        polyorder = params['degree']
    ) 
    return aggregated
    
def label_svs_with_null_statistics(svs, aggregated_over_length_quantiles):     
    sv_length_quantile = np.array(aggregated_over_length_quantiles['sv length quantile'])
    estimated_singleton_probability_under_null = np.array(aggregated_over_length_quantiles['estimated singleton probability under null'])
    mean_isSingleton_null = estimated_singleton_probability_under_null
    variance_isSingleton_null = estimated_singleton_probability_under_null*(1 - estimated_singleton_probability_under_null)    
    null_statistics = pd.DataFrame(data={
        'sv length quantile': sv_length_quantile, 
        'mean(isSingleton) under null': mean_isSingleton_null,
        'variance(isSingleton) under null': variance_isSingleton_null
    })
    return pd.merge(svs, null_statistics)
    
def label_svs_with_score_quantiles(df, score, operator, params): 
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")        

        array_of_quantiles = params['array_of_quantiles']
        starts = array_of_quantiles[:-1]
        ends = array_of_quantiles[1:]
        quantile_labels = [f'{start} - {end}' for start, end in zip(starts, ends)]

        df[f'({operator} {score}) quantile'], bins = pd.qcut(
            df[f'{operator} {score}'],
            q = array_of_quantiles, 
            labels = quantile_labels, 
            retbins = True,
    #         duplicates='drop'
        )

        return df

def aggregate_over_score_quantiles(df, score, operator): 
    groups = df.groupby([f'({operator} {score}) quantile'])
    aggregated = groups.agg({
        'sv is singleton': ['sum', 'count'],
        'mean(isSingleton) under null': ['sum', 'count'],
        'variance(isSingleton) under null': ['sum', 'count'],
        f'{operator} {score}': ['mean', 'std']
    })
    aggregated = aggregated.reset_index()
    aggregated['observed F'] = aggregated[('sv is singleton', 'sum')]/aggregated[('sv is singleton', 'count')]
    aggregated['mean(F) under null'] = aggregated[('mean(isSingleton) under null', 'sum')]/aggregated[('mean(isSingleton) under null', 'count')]
    aggregated['variance(F) under null'] = aggregated[('variance(isSingleton) under null', 'sum')]/(aggregated[('variance(isSingleton) under null', 'count')]**2)
    aggregated['observed F, relative to null'] = aggregated['observed F']/aggregated['mean(F) under null']
    aggregated['zscore(F)'] = (aggregated['observed F'] - aggregated['mean(F) under null'])/np.sqrt(aggregated['variance(F) under null'])
    return aggregated

def resample_aps(source, score, operator): 
  svs = get_svs(source, do_resample=True)
  params = parameters_to_estimate_singleton_probability_under_null[source]  
  svs = label_svs_with_length_quantiles(svs, params)
  aggregated_over_length_quantiles = aggregate_over_length_quantiles(svs, params)
  svs = label_svs_with_null_statistics(svs, aggregated_over_length_quantiles)
  svs = label_svs_with_score_quantiles(svs, score, operator, params)
  aggregated_over_score_quantiles = aggregate_over_score_quantiles(svs, score, operator)
  aps = aggregated_over_score_quantiles[[
    f'({operator} {score}) quantile',
    f'{operator} {score}',
    'observed F, relative to null'
  ]]
  aps.columns = [' '.join(col[::-1]).strip() for col in aps.columns.values]
  return aps

def ci_lower(xs):
  xs = np.array(xs)
  return np.percentile(xs, 2.5) 

def ci_upper(xs): 
  xs = np.array(xs)
  return np.percentile(xs, 97.5) 

def bootstrap_aps(source, score, number_bootstrap_samples, operator='min'): 
  aps_samples = []
  for i in tqdm(range(number_bootstrap_samples), desc=f'{source}; {score}'):
    aps_samples.append(resample_aps(source, score, operator))
  aps_samples = pd.concat(aps_samples, axis=0, ignore_index=True)

  groups = aps_samples.groupby([f'({operator} {score}) quantile'])
  aggregated = groups.agg({
      'observed F, relative to null': [
        'mean', 
        ('number of bootstrap samples', 'count'),
        ('ci_lower', ci_lower), 
        ('ci_upper', ci_upper), 
      ],
  })
  aggregated = aggregated.reset_index()
  return aggregated

def save_aps(number_bootstrap_samples): 
  sources = ['gnomAD', 'CCDG', '1000G']
  scores = ['negative new chen zscore']
  for kmer_size in [3, 5, 7]: 
      for train_set_label in ['coding', 'noncoding', 'chenWindows']: 
          scores.append(f'N_bar_{kmer_size}_{train_set_label}')
  aps = defaultdict(lambda: {})
  for source in sources: 
    for score in scores: 
      aps[source][score] = bootstrap_aps(source, score, number_bootstrap_samples)
  aps = dict(aps)
  with open('aps.pkl', 'wb') as handle:
      pickle.dump(aps, handle, protocol=pickle.HIGHEST_PROTOCOL)

if __name__ == '__main__':
  save_aps(number_bootstrap_samples=int(sys.argv[1]))