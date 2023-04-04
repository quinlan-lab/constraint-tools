# pysam API: 
# https://github.com/pysam-developers/pysam/blob/b82cbcae22c088e64fdb58f8acaf1e9773c7b088/pysam/libctabix.pyx
import pysam

import color_traceback
from compute_Nbar_Nobserved import compute_Nbar_Nobserved 

def compute_zscores_on_window(window, model, log=True):
  with pysam.TabixFile(model['mutations']) as mutations, pysam.FastaFile(model['genome']) as genome:
    N_bar, N_observed, N_mean_null, N_variance_null = compute_Nbar_Nobserved(window, model, mutations, genome, log, verbose=True) 
  return N_bar, N_observed, N_mean_null, N_variance_null

