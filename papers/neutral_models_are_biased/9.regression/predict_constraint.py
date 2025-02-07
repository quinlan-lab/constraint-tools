from sklearn.metrics import precision_recall_curve
import seaborn as sns
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import numpy as np

from generate_data import plot_xs

NUMBER_EXAMPLES_MIN = 100

def filter_by_bin(df, bin_center, bin_width):
    left = bin_center - bin_width / 2
    right = bin_center + bin_width / 2
    df = df[
        (df['x'] >= left) & 
        (df['x'] < right)
    ]
    return df 

def normal_distribution(x, mu, sigma):
    return (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

def plot_residual_distributions(ax, df, standardized, xlim, model_type, title, ylabel, ylim):
    xlabel = f'standardized_residuals_{model_type}Model' if standardized else f'residuals_{model_type}Model'    
    xlabel_plot = f'standardized residual ({model_type} model)' if standardized else f'residual ({model_type} model)'    

    bins = np.linspace(xlim[0], xlim[1], 100)         
    bin_centers = (bins[1:] + bins[:-1]) / 2
    bin_width = bins[1] - bins[0]

    palette = {False: 'blue', True: 'orange'}
    alpha = 0.5  

    sns.histplot(
        data=df, 
        x=xlabel, 
        hue='constrained', 
        kde=False, 
        element='step', 
        bins=bins, 
        ax=ax, 
        palette=palette, 
        alpha=alpha
    )

    # Plot normal distribution
    x = bin_centers
    if standardized:
        y = normal_distribution(x, mu=0, sigma=1)
    else:
        # assumes the x in df are in the vicinity of x0, 
        # so that the average y value of the negative examples is a good estimate of TRUE_RATE(x0)
        lambda_ = df[df['constrained'] == False]['y'].mean()

        y = normal_distribution(x, mu=0, sigma=np.sqrt(lambda_)) 
    negative_class_count = len(df[df['constrained'] == False])
    y_scaled = y * negative_class_count * bin_width  
    normal_line, = ax.plot(x, y_scaled, color='red', lw=2)

    ax.set_xlabel(xlabel_plot)
    ax.set_ylabel(ylabel)
    ax.set_yscale('log')  
    ax.set_xlim(xlim) 
    ax.set_ylim(ylim)
    ax.set_title(title)

    # Manually create legend handles and labels
    handles = [
        Patch(color=palette[False], alpha=alpha),
        Patch(color=palette[True], alpha=alpha),
        normal_line
    ]
    labels = [
        'neutral', 
        'constrained', 
        'standard normal dist.' if standardized else 'normal dist.'
    ]
    ax.legend(handles=handles, labels=labels, prop={'size': 15}, loc='lower center')

def plot_pr_curve(ax, df, model_type, standardized):
    # using standardized_residuals appear to yield larger auPRC than raw residuals
    # c.f., section entitled "Model bias is responsible for poor genome-wide performance" in this notebook
    if standardized:
        precision, recall, _ = precision_recall_curve(df['constrained'], df[f'standardized_residuals_{model_type}Model'])
    else:
        precision, recall, _ = precision_recall_curve(df['constrained'], df[f'residuals_{model_type}Model'])
    num_examples = len(df)

    if num_examples > NUMBER_EXAMPLES_MIN:
        ax.plot(recall, precision, label=f'{model_type} model')

def compute_delta_sigma(df):
    df_neutral = df[df['constrained'] == False]
    df_constrained = df[df['constrained'] == True]
    delta = df_neutral['y'].mean() - df_constrained['y'].mean()
    sigma = df_neutral['y'].std()
    return delta, sigma

def annotate_plot_residual_distributions(ax, delta, sigma): 
    ax.axvline(x=0, color='green', linestyle='--', lw=2)
    ax.axvline(x=delta, color='green', linestyle='--', lw=2)
    ax.annotate(
        text='', 
        xy=(0, ax.get_ylim()[1] / 2), 
        xytext=(delta, ax.get_ylim()[1] / 2), 
        arrowprops=dict(arrowstyle='<->', color='green', lw=2)
    )
    ax.text(
        x=delta / 2, 
        y=(ax.get_ylim()[1] / 2) * 1.5, 
        s=r'$\Delta$', 
        ha='center', 
        va='center', 
        color='green', 
        fontsize=15
    )

    ax.annotate(
        text='', 
        xy=(-sigma, ax.get_ylim()[1] / 4), 
        xytext=(sigma, ax.get_ylim()[1] / 4), 
        arrowprops=dict(arrowstyle='<->', color='red', lw=2)
    )
    ax.text(
        x=0, 
        y=(ax.get_ylim()[1] / 4) * 1.5, 
        s=r'$2\sigma$', 
        ha='center', 
        va='center', 
        color='red', 
        fontsize=15
    )
    
def plot_pr_curve_wrapper(df, model_types, positive_fraction, xlim_residual, ylim_residual, standardized, number_examples, bin_widths): 
    bin_centers = [0, 1, 2]

    # duplicate the first element so we can plot the distribution of residuals for the smallest bin width:
    bin_widths.insert(0, bin_widths[0]) 
    
    number_rows = len(bin_centers)
    number_columns = len(bin_widths) + 1 # +1 for an additional column that shows the distribution of feature values
    fig, axes = plt.subplots(number_rows, number_columns, figsize=(5*number_columns, 5*number_rows))

    for i, bin_center in enumerate(bin_centers):
        plot_xs(axes[i, 0], number_examples, xlim=(-0.5, 2.5), yscale='linear')
        for j, bin_width in enumerate(bin_widths):
            df_filtered = filter_by_bin(df, bin_center, bin_width)

            j += 1 # skip the first column
            ax = axes[i, j]
            if j == 1:
                if bin_width < 0.1 + 1e-6: 
                    delta, sigma = compute_delta_sigma(df_filtered) 
                else: 
                    raise ValueError('bin_width too large to reliably compute delta (depletion of SNVs) and sigma (std of predicted SNV counts)')
                SNR = 1 + (delta/sigma)**2 # signal-to-noise ratio
                plot_residual_distributions(
                    ax,
                    df_filtered, 
                    standardized=standardized, 
                    xlim=xlim_residual,
                    model_type='quadratic', 
                    title=(
                        f'SNR = {SNR:.2f}\n' 
                        f'feature-bin width = {bin_width}'
                    ),
                    ylabel=(
                        'number of windows\n'
                        f'(feature-bin center = {bin_center})'
                    ), 
                    ylim=ylim_residual
                )
                if not standardized:
                    annotate_plot_residual_distributions(ax, delta, sigma)

            else: 
                axes[i, 0].axvspan(bin_center - bin_width / 2, bin_center + bin_width / 2, color='red', alpha=0.3)

                for model_type in model_types:                
                    plot_pr_curve(ax, df_filtered, model_type, standardized)
                ax.plot([0, 1], [positive_fraction, positive_fraction], color='black', lw=2, linestyle=':', label='random classifier')
                ax.set_xlabel('Recall')
                ax.set_ylabel('Precision')
                ax.set_title(f'feature-bin width: {bin_width}')
                ax.legend(prop={'size': 15})  
                ax.set_xlim(0, 1)
                ax.set_ylim(0, 1)
                ax.grid(True)
                
    plt.tight_layout()
    plt.show()

def plot_residual_distributions_all_models(df, standardized, xlim, model_types): 
    fig, axes = plt.subplots(1, len(model_types), figsize=(7*len(model_types), 5))

    for i, model_type in enumerate(model_types):
        ax = axes[i]
        plot_residual_distributions(
            ax, 
            df, 
            standardized=standardized, 
            xlim=xlim, 
            model_type=model_type, 
            title=f'{model_type} model', 
            ylabel='number of examples', 
            ylim=(1, 1e5)
        )

    plt.tight_layout()
    plt.show()

