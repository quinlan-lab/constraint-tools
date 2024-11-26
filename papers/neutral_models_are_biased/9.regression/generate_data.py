import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def generate_xs(number_examples):
    # Sample from a univariate normal distribution
    MEAN = 0
    STD = 1
    xs = np.random.normal(loc=MEAN, scale=STD, size=number_examples)
    return xs

def plot_xs(number_examples):
    xs = generate_xs(number_examples)

    plt.figure(figsize=(6, 6))
    plt.hist(xs, bins=100)
    plt.xlabel('genomic feature')
    plt.ylabel('number of windows')
    plt.yscale('log')
    # plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: '{:.0e}'.format(y)))
    plt.show()

def compute_true_params(A, B, C):
    return { 
        'keys': ['A', 'B', 'C'], 
        'values': [A, B, C]
    }

def compute_true_rate(true_params):
    A, B, C = true_params['values'] 
    true_rate = lambda x: np.exp(A + B*x + C*x**2) # rate function
    return true_rate

def generate_ys(true_params, xs):
    true_rate = compute_true_rate(true_params)
    ys = np.random.poisson(lam=true_rate(xs))
    return ys 

def plot_ys(true_params, number_examples):
    xs = generate_xs(number_examples)
    ys = generate_ys(true_params, xs)

    plt.figure(figsize=(12, 6))
    plt.hist(ys, bins=1000)
    plt.xlabel('SNV counts (before selection)')
    plt.ylabel('number of windows')
    plt.xlim(0, 400)

def compute_y_pos_1(xs, ys, num_pos):
    y_depletion = 15 # informed by experiment
    y_pos = np.maximum(0, ys[:num_pos] - y_depletion)
    return y_pos

# mimic the scenario where enhancers, say, are enriched in GC content: 
def compute_y_pos_2(xs, ys, num_pos): 
    fractional_reduction_in_y = 0.2
    x_factor = (xs - xs.min()) / (xs.max() - xs.min())
    y_pos = ys[:num_pos] * (1 - fractional_reduction_in_y * x_factor[:num_pos])
    return y_pos

def generate_xs_ys_with_selection(true_params, number_examples, positive_fraction, compute_y_pos):
    xs = generate_xs(number_examples)
    ys = generate_ys(true_params, xs)

    # Declare a random collection of the examples to be positive examples, 
    # and reduce their y values, to mimic the effect of negative selection 
    num_pos = int(positive_fraction*number_examples)
    y_pos = compute_y_pos(xs, ys, num_pos)

    # Declare the remaining examples to be neutral (negative) examples, and do not change their y values
    y_neg = ys[num_pos:]

    ys = np.concatenate([y_pos, y_neg])
    constrained = num_pos*[True] + (number_examples-num_pos)*[False]

    data = pd.DataFrame({'x': xs, 'y': ys, 'constrained': constrained})

    import seaborn as sns

    plt.figure(figsize=(12, 6))
    sns.scatterplot(data=data, x='x', y='y', hue='constrained', palette='viridis', alpha=0.5)
    plt.xlabel('genomic feature')
    plt.ylabel('SNV counts')
    plt.yscale('log')

    return data
