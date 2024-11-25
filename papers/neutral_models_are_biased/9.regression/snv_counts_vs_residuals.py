import pandas as pd
import matplotlib.pyplot as plt

def plot_snv_counts_vs_residuals(df, model_type):
    df = pd.DataFrame({ 
        'y': df['y'],
        f'predicted_y_{model_type}Model': df[f'predicted_y_{model_type}Model'],
        f'standardized_residuals_{model_type}Model_bin': pd.qcut(
            df[f'standardized_residuals_{model_type}Model'], 25, labels=None
        ), 
    })

    df_grouped = ( 
        df
        .groupby(f'standardized_residuals_{model_type}Model_bin')
        .agg({
            'y': 'mean',
            f'predicted_y_{model_type}Model': 'mean',
        })
    )

    df_grouped['y_all_bins'] = df['y'].mean()

    df_grouped.plot(
        y=['y', f'predicted_y_{model_type}Model', 'y_all_bins'], 
        kind='line', 
        figsize=(8, 5)
    )

    plt.xlabel(f'standardized residuals ({model_type}Model)')
    plt.ylabel('mean value')
    plt.legend()
    plt.xticks(rotation=45)
    plt.ylim(100, 300)
    plt.show()

def plot_snv_counts_vs_residuals_all_models(df, model_types):
    for model_type in model_types:
        plot_snv_counts_vs_residuals(df, model_type)

