from singleton import update_types
import json 

def read_model(filename): 
  with open(filename, 'r') as fh:
    d = json.load(fh) 
    if ('singletonCounts' in d.keys()) and ('singletonProbabilities' in d.keys()): 
      d['singletonCounts'] = update_types(d['singletonCounts'])
      d['singletonProbabilities'] = update_types(d['singletonProbabilities'])
    return d

def test_read_model(): 
  import sys 
  from colorize import print_json, print_string_as_info, print_string_as_info_dim
  import numpy as np

  model_filename = sys.argv[1]

  with open(model_filename, 'r') as fh:
    d = json.load(fh) 
    print_string_as_info('first item in singletonCounts prior to updating types: ')
    for k, v in d['singletonCounts'].items(): 
      assert isinstance(k, str)
      assert isinstance(v, list)
      print(k, type(k))
      print(v, type(v))
      break
    print_string_as_info('first item in singletonProbabilities prior to updating types: ')
    for k, v in d['singletonProbabilities'].items(): 
      assert isinstance(k, str)     
      assert isinstance(v, list)     
      print(k, type(k))
      print(v, type(v))
      break

  print_string_as_info('updating types of model....')
  model = read_model(model_filename) 

  print_string_as_info('mutations:')
  print_string_as_info_dim(model['mutations'])
  print(type(model['mutations']))

  print_string_as_info('kmerSize:')
  print(model['kmerSize'])
  print(type(model['kmerSize']))

  print_string_as_info("kmerCounts['AAAAAAA']:")
  print(model['kmerCounts']['AAAAAAA'])
  print(type(model['kmerCounts']['AAAAAAA']))

  print_string_as_info("singletonCounts:")
  for k, v in model['singletonCounts'].items(): 
    assert isinstance(k, int)     
    assert isinstance(v, np.ndarray)     
    print(k, type(k))    
    print(v, type(v))
    break

  print_string_as_info("singletonProbabilities:")
  for k, v in model['singletonProbabilities'].items(): 
    assert isinstance(k, int)     
    assert isinstance(v, np.ndarray)     
    print(k, type(k))
    print(v, type(v))
    break

if __name__ == '__main__': 
  test_read_model() 

