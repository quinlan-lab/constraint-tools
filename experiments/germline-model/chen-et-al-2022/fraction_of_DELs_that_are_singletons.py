import pandas as pd 
import numpy as np 
from scipy.signal import savgol_filter
import warnings 

parameters_to_estimate_singleton_probability_under_null = { 
    'gnomAD': dict(
      number_length_quantiles=200, 
      savgol_window_length=51, 
      degree=3,
      array_of_quantiles=[0.000, 0.02, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.98, 0.99, 1.00]
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
    filename = f'{CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/{source}-noncoding-svs-chen-mchale.bed'
    df = pd.read_csv(filename, sep='\t')
    df['negative chen zscore'] = -df['chen_zscore']
    
    df['K_bar'] = pd.to_numeric(df['K_bar'], errors='coerce') # convert '.' to 'NaN'
#     df = df.dropna(subset=['K_bar']) # drop windows for which K_bar == NaN
    df['negative K_bar'] = -df['K_bar']
    
    df = df[['sv_id', 'sv_length', 'alt_allele_count', 'N_bar', 'negative K_bar', 'negative chen zscore']]
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
    aggregated = groups.agg({
        'N_bar': ['min', 'mean'],
        'negative K_bar': ['min', 'mean'],
        'negative chen zscore': ['min', 'mean']
    })
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
        'variance(isSingleton) under null': ['sum', 'count']
    })
    aggregated = aggregated.reset_index()
    aggregated['observed F'] = aggregated[('sv is singleton', 'sum')]/aggregated[('sv is singleton', 'count')]
    aggregated['mean(F) under null'] = aggregated[('mean(isSingleton) under null', 'sum')]/aggregated[('mean(isSingleton) under null', 'count')]
    aggregated['variance(F) under null'] = aggregated[('variance(isSingleton) under null', 'sum')]/(aggregated[('variance(isSingleton) under null', 'count')]**2)
    aggregated['Adjusted Proportion of Singletons'] = aggregated['observed F'] - aggregated['mean(F) under null']
    aggregated['zscore(F)'] = (aggregated['observed F'] - aggregated['mean(F) under null'])/np.sqrt(aggregated['variance(F) under null'])
    return aggregated

