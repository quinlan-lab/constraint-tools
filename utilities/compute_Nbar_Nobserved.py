import numpy as np 

from snvs import fetch_SNVs
from get_p0s_p1s_p2s_p3s import get_p0s_p1s_p2s_p3s

def get_N_observed(window, genome, mutations, model):
  return len(fetch_SNVs(mutations, genome, window['region'], meta=model))

# https://github.com/quinlan-lab/constraint-tools/blob/main/define-model/germline-model.ipynb
def get_N_mean_variance_null(window, genome, model, log): 
  p0s_p1s_p2s_p3s = get_p0s_p1s_p2s_p3s(window, genome, model, log)

  mean_Ns = np.array([0*p0 + 1*p1 + 2*p2 + 3*p3 for p0, p1, p2, p3 in p0s_p1s_p2s_p3s])
  mean_N = np.sum(mean_Ns)

  mean_N2s = np.array([(0**2)*p0 + (1**2)*p1 + (2**2)*p2 + (3**2)*p3 for p0, p1, p2, p3 in p0s_p1s_p2s_p3s])
  variance_N = np.sum(mean_N2s - np.square(mean_Ns))

  return mean_N, variance_N

def compute_Nbar_Nobserved(window, model, mutations, genome, log):
  N_mean_null, N_variance_null = get_N_mean_variance_null(window, genome, model, log)
  N_observed = get_N_observed(window, genome, mutations, model)
  N_bar = (N_observed - N_mean_null)/np.sqrt(N_variance_null)
  return N_bar, N_observed

