from collections import defaultdict 
import numpy as np 

from colorize import print_unbuffered, print_string_as_info

# https://stackoverflow.com/a/43268377/6674256
# https://docs.python.org/3/library/collections.html#collections.defaultdict.__missing__
class Singleton_Counts(defaultdict):
  def __missing__(self, SNV_count):
    self[SNV_count] = singleton_counts = np.array([0]*(SNV_count+1))
    return singleton_counts

def add_singleton_counts(x, y):
  for SNV_count in y.keys():
    x[SNV_count] += y[SNV_count]
  return x

def defaultdict_to_dict(singleton_counts): 
  return {
    k: v.tolist() # numpy array to list 
    for k, v in singleton_counts.items()
  }

def test_data_structure(): 
  print_unbuffered('')
  print_string_as_info('********* testing data structure ****************')
  print_unbuffered('')

  singleton_counts = Singleton_Counts()
  print(singleton_counts) 
  singleton_counts[4][3] += 1 
  print('after addition of SNV count and singleton count:') 
  print(singleton_counts) 

def test_combine_singleton_counts(): 
  print_unbuffered('')
  print_string_as_info('********* testing function to combine two data structures ****************')
  print_unbuffered('')

  x = Singleton_Counts()
  x[5][3] += 1
  x[4][1] += 1
  x[5][1] += 1
  print('x:', x)

  y = Singleton_Counts()
  y[5][3] += 1
  y[6][2] += 1
  y[1][1] += 1
  print('y:', y)

  print('combine x and y:') 
  print(add_singleton_counts(x, y)) 

  print('x:', x)

if __name__ == '__main__': 
  test_data_structure()
  test_combine_singleton_counts()


