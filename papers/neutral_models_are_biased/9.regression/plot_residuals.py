import pandas as pd
import matplotlib.pyplot as plt

def plot_residuals_vs_feature(ax, df, standardized, model_type):
    df = df[df['constrained'] == False]

    df['x_bin_midpoints'] = (pd
        .cut(df['x'], bins=20, labels=None)
        .apply(lambda x: (x.right + x.left) / 2)
    )

    residual_label = f'standardized_residuals_{model_type}Model' if standardized else f'residuals_{model_type}Model'    
    residual_label_plot = f'standardized residuals' if standardized else f'residuals'    

    average_residual = df.groupby('x_bin_midpoints')[residual_label].mean()  

    ax.scatter(df['x'], df[residual_label], alpha=0.7, label=residual_label_plot)
    ax.plot(average_residual.index, average_residual.values, color='red', lw=2, label=f'average') 
    ax.axhline(0, color='black', linestyle='--')

    ax.set_title(f'{model_type} model')
    ax.set_xlabel('genomic feature')
    ylim = (-10, 10) if standardized else (-100, 100) 
    ax.set_ylim(ylim)
    yticks = [-10, -5, 0, 5, 10] if standardized else [-100, -50, 0, 50, 100] 
    ax.set_yticks(yticks)
    ax.set_xlim(-5, 5)
    ax.legend(prop={'size': 20})

def plot_residuals_vs_feature_all_models(df, standardized, model_types): 
    fig, axes = plt.subplots(1, len(model_types), figsize=(7*len(model_types), 5))

    for i, model_type in enumerate(model_types):
        ax = axes[i]
        plot_residuals_vs_feature(
            ax, 
            df, 
            standardized, 
            model_type, 
        )

    plt.tight_layout()
    plt.show()
