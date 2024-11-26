import pandas as pd
import matplotlib.pyplot as plt

def plot_snv_counts_vs_residuals(df, model_type):
    df = pd.DataFrame({ 
        'observed SNV counts': df['y'],
        f'predicted SNV counts ({model_type} model)': df[f'predicted_y_{model_type}Model'],
        f'standardized residuals bin ({model_type} model)': pd.qcut(
            df[f'standardized_residuals_{model_type}Model'], 25, labels=None
        ), 
    })

    df_grouped = ( 
        df
        .groupby(f'standardized residuals bin ({model_type} model)')
        .agg({
            'observed SNV counts': 'mean',
            f'predicted SNV counts ({model_type} model)': 'mean',
        })
    )

    df_grouped['SNV counts (all bins)'] = df['observed SNV counts'].mean()

    df_grouped.plot(
        y=['observed SNV counts', f'predicted SNV counts ({model_type} model)', 'SNV counts (all bins)'], 
        kind='line', 
        figsize=(8, 5)
    )

    plt.xlabel(f'standardized residuals ({model_type} model)')
    plt.ylabel('mean value')
    plt.legend()
    plt.xticks(rotation=45)
    plt.ylim(100, 300)
    plt.show()

def plot_snv_counts_vs_residuals_all_models(df, model_types):
    for model_type in model_types:
        plot_snv_counts_vs_residuals(df, model_type)

