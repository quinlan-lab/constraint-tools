import numpy as np 

def ravel(sequences):  
  number_alleles, number_sites = sequences.shape
  return np.reshape(sequences, number_alleles*number_sites, order='C')

def unravel(embeddings, number_alleles, number_sites): 
  dim, embedding_dimension = embeddings.shape
  assert dim == number_alleles*number_sites
  return np.reshape(embeddings, (number_alleles, number_sites, embedding_dimension), order='C')   
