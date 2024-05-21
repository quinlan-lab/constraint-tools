import matplotlib.pyplot as plt

def length_to_string(length): 
    if length < 1000: 
        return f'{length}bp'
    elif length < 1000000: 
        return f'{int(length/1000)}kb'
    else: 
        return f'{int(length/1000000)}Mb'
    
def compute_limits(df, feature, mean_factor, std_factor):
    mean = df[feature].mean()
    scaled_mean = mean_factor*mean
    std = df[feature].std()
    scaled_std = std_factor*std
    return scaled_mean - scaled_std, scaled_mean + scaled_std

def slice_feature_space(df, conditional_features_and_lims): 
  for conditional_feature, lim in conditional_features_and_lims: 
    print(f'conditioning on {conditional_feature} in [{lim[0]}, {lim[1]}]')
    df = df[
      (df[conditional_feature] > lim[0]) & 
      (df[conditional_feature] < lim[1])
    ]
  return df 

def plot_feature_distribution(df, feature, xlabel, lim=None): 
    x = df[feature] 

    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    ax.hist(
        x, 
        density=True, 
        bins=100,
        histtype='stepfilled', 
        alpha=0.2, 
        label='observed', 
        color='black'
    )

    if lim is not None:
        ax.axvspan(lim[0], lim[1], alpha=0.2, color='red', label='conditioned values')

    ax.set_xlabel(xlabel)
    ax.set_ylabel('probability density')
    ax.set_yscale('linear')
    ax.set_xlim(x.min(), x.max())
    ax.legend()
    
