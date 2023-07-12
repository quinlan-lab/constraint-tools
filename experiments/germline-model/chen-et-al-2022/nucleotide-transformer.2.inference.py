CONSTRAINT_TOOLS = '/scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools'
CONSTRAINT_TOOLS_DATA = '/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools'

import pandas as pd 

pd.set_option('display.max_columns', 50)

filename_prefix = f'{CONSTRAINT_TOOLS_DATA}/benchmark-genome-wide-predictions/chen-et-al-2022/enhancer-characteristics-enrichment-subset' 

def get_windows_scores_annotations():
  filename = f'{filename_prefix}.bed'
  df = pd.read_csv(filename, sep='\t')
  df = df[df['window overlaps merged_exon'] == False]
  return df

windows_scores_annotations_noncoding = get_windows_scores_annotations()

##########################################################################

# based upon:
# https://github.com/instadeepai/nucleotide-transformer/blob/main/examples/inference.ipynb

# I installed nucleotide-transformer using: 
# pip install git+https://github.com/instadeepai/nucleotide-transformer@main

# Another way to use nucleotide-transformer: 
# https://huggingface.co/InstaDeepAI/nucleotide-transformer-500m-human-ref?text=ACCTGA%3Cmask%3EAACTGAGTC

import haiku as hk
import jax
import jax.numpy as jnp
from nucleotide_transformer.pretrained import get_pretrained_model

model_name = '500M_human_ref' #@param['500M_human_ref', '500M_1000G', '2B5_1000G', '2B5_multi_species']

# TODO: change these parameters: 
# with the old values for the following parameters, 
# one site takes 20secs on CPU => 1000 sites takes 5hrs on CPU 
SEQUENCE_LENGTH = 11 # 501 # 5977 # must be odd
NUMBER_TOKENS_PER_SEQUENCE = 15 # 100 # 1000 # split sequence into book-ended 6-mer tokens and add CLS token

parameters, forward_fn, tokenizer, config = get_pretrained_model(
    model_name=model_name,
    mixed_precision=False,
    embeddings_layers_to_save=(20,), # the layers at which you'd like to collect embeddings (e.g. (5, 10, 20) to get embeddings at layers 5, 10 and 20)
    attention_maps_to_save=((1, 4), (7, 18)), # the attention maps you´d like to collect (e.g. ((1,4), (7,18)) to get attention maps corresponding to layer 1 head number 4 and layer 7 head number 18)
    max_positions=NUMBER_TOKENS_PER_SEQUENCE # we recommend keeping this number as small as possible for optimized memory and inference time
)
forward_fn = hk.transform(forward_fn)

# L91 of $HOME/.conda/envs/constraint-tools/lib/python3.9/site-packages/nucleotide_transformer/model.py : 
assert NUMBER_TOKENS_PER_SEQUENCE <= config.max_positions 

# for inference: 
random_key = jax.random.PRNGKey(0)

##########################################################################

import numpy as np 

import sys
sys.path.append(f'{CONSTRAINT_TOOLS}/utilities')

from kmer import middle_index

def set_center_nucleotide(sequence, new_center_nucleotide): 
  sequence = list(sequence)
  sequence[middle_index(sequence)] = new_center_nucleotide
  return ''.join(sequence)
   
def test_set_center_nucleotide(): 
   sequence = 'AAGCT'
   new_center_nucleotide = 'T'
   new_sequence = f'AA{new_center_nucleotide}CT'
   assert set_center_nucleotide(sequence, new_center_nucleotide) == new_sequence

test_set_center_nucleotide() 

##########################################################################

import pysam 

from read_model import read_model
from kmer import fetch_kmers
from pack_unpack import pack 
from ravel_unravel import ravel, unravel
from bases import BASES 

MCHALE_MODEL_FILENAME = f"{CONSTRAINT_TOOLS}/dist/model-germline-grch38-Nonly.kmerSize-3.trainSet-noncoding.json"
MCHALE_MODEL = read_model(MCHALE_MODEL_FILENAME)

def do_inference_on_window(window): 
  chromosome = window['chromosome']
  start = window['start']
  end = window['start'] + 3 # TODO: replace with window['end']
  region = pack(chromosome, start, end)
  
  sequences = []
  with pysam.FastaFile(MCHALE_MODEL['genome']) as genome:
    for sequence in fetch_kmers(region, genome, kmer_size=SEQUENCE_LENGTH, log=True): # upper case
      alleles = []
      for new_center_nucleotide in BASES:
        allele = set_center_nucleotide(sequence, new_center_nucleotide)
        alleles.append(allele)
      sequences.append(alleles)
  sequences = np.array(sequences)
  sequences = sequences.T # [A, C, G, T] X [site1, site2, ...]
  sequences_raveled = ravel(sequences)

  # Tokenize sequence(s)
  # The `tokens_str` variable shows how sequence(s) have been split into tokens. 
  # The token list will be padded to size `max_positions`.
  tokens_ids = [b[1] for b in tokenizer.batch_tokenize(list(sequences_raveled))] # batch_tokenize expects "List[str]"
  # tokens_str = [b[0] for b in tokenizer.batch_tokenize(list(sequences_raveled))] # batch_tokenize expects "List[str]"
  tokens = jnp.asarray(tokens_ids, dtype=jnp.int32)

  # Perform inference
  # The first time you query this cell, it will be slower than usual because of the computation graph compilation.        
  outs = forward_fn.apply(parameters, random_key, tokens)  

  # Retrieve embeddings
  embeddings = outs["embeddings_20"][:, 1:, :]  # removing CLS token
  padding_mask = jnp.expand_dims(tokens[:, 1:] != tokenizer.pad_token_id, axis=-1)
  masked_embeddings = embeddings * padding_mask  
  sequences_lengths = jnp.sum(padding_mask, axis=1)
  mean_embeddings = jnp.sum(masked_embeddings, axis=1) / sequences_lengths

  mean_embeddings = unravel(mean_embeddings, number_alleles=len(BASES), number_sites=end-start)
  mean_embeddings = np.array((mean_embeddings))

  return mean_embeddings

def save_inference(df): 
  for i, window in df.iterrows(): 
    window.to_pickle(f'{filename_prefix}.{i}.pkl')  
    mean_embeddings = do_inference_on_window(window)
    np.save(f'{filename_prefix}.{i}.npy', mean_embeddings)
    print('saved windows and embeddings to:')
    print(f'{filename_prefix}.{i}.*')

save_inference(windows_scores_annotations_noncoding)