from kmer import (
  fetch_kmers,
  compute_possible_ALT_states_core
)

def get_substitution_probability(model, kmer, ALT_multiplicity): 
  return sum(
    model['kmerProbabilities'][kmer][alt_state] 
    for alt_state in compute_possible_ALT_states_core(kmer, ALT_multiplicity)
  )

# https://github.com/quinlan-lab/constraint-tools/blob/main/define-model/germline-model.ipynb
def get_p0s_p1s_p2s_p3s(window, genome, model, log=True):
  p0s_p1s_p2s_p3s = []
  for kmer in fetch_kmers(window['region'], genome, model['kmerSize'], log): 
    p1_p2_p3 = [get_substitution_probability(model, kmer, ALT_multiplicity) for ALT_multiplicity in [1, 2, 3]]
    p0 = 1 - sum(p1_p2_p3) 
    p0s_p1s_p2s_p3s.append([p0] + p1_p2_p3)
  return p0s_p1s_p2s_p3s

