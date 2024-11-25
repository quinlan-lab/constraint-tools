import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import statsmodels.api as sm

from generate_data import compute_true_rate

# Mute SettingWithCopyWarning
pd.options.mode.chained_assignment = None

def plot_y(df, model_type, true_params):
    df = df[df['constrained'] == False]

    x_ = np.linspace(df['x'].min(), df['x'].max(), 100)
    true_rate = compute_true_rate(true_params)
    y_ = true_rate(x_)

    plt.figure(figsize=(12, 6))
    plt.plot(df['x'], df['y'], 'o', alpha=0.5, label='y')
    plt.plot(df['x'], df[f'predicted_y_{model_type}Model'], 'o', label=f'Learned rate under {model_type} model')
    plt.plot(x_, y_, label='True rate', color='black')
    plt.yscale('log')
    plt.xlabel('x')
    plt.legend(prop={'size': 25})
    plt.xlim(-5, 5)
    plt.ylim(1e1, 1e4)
    plt.show()

def fit_poisson_model(df, model_type, true_params):
    df_neg = df[df['constrained'] == False]

    if model_type == 'constant':
        x_model_neg = np.ones((df_neg.shape[0], 1))
        x_model = np.ones((df.shape[0], 1))  
    elif model_type == 'linear':
        x_model_neg = sm.add_constant(df_neg['x'])
        x_model = sm.add_constant(df['x'])
    elif model_type == 'quadratic':
        df_neg['x2'] = df_neg['x']**2
        x_model_neg = sm.add_constant(df_neg[['x', 'x2']])
        df['x2'] = df['x']**2
        x_model = sm.add_constant(df[['x', 'x2']])
    else:
        raise ValueError(f"Unknown model type: {model_type}")
            
    neutral_model = sm.Poisson(df_neg['y'], x_model_neg).fit(disp=False)

    print(f'{model_type} model of lambda:')
    estimated_params_keys = ['alpha', 'beta', 'gamma']
    estimated_params_values = neutral_model.params.values
    for i, estimated_param_value in enumerate(estimated_params_values):
        print(f'{estimated_params_keys[i]}: {estimated_param_value:.2f} ({true_params["keys"][i]}: {true_params["values"][i]})')
    print()
    
    df[f'predicted_y_{model_type}Model'] = neutral_model.predict(x_model)

    df[f'residuals_{model_type}Model'] = df[f'predicted_y_{model_type}Model'] - df['y']
    df[f'standardized_residuals_{model_type}Model'] = df[f'residuals_{model_type}Model']/np.sqrt(df[f'predicted_y_{model_type}Model'])
        
    plot_y(df, model_type, true_params)

    return df

def fit_poisson_model_wrapper(df, model_types, true_params):
    for model_type in model_types:
        df = fit_poisson_model(df, model_type, true_params)
    return df
