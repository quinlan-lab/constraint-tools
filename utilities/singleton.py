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

def make_serializable(singleton_data): 
  return {
    k: v.tolist() # numpy array to list 
    for k, v in singleton_data.items()
  }

def dict_to_defaultdict(singleton_data):
  singleton_data = {
    k: np.array(v) # list to numpy array 
    for k, v in singleton_data.items()
  }
  return Singleton_Counts(None, singleton_data)

def test_dict_to_defaultdict(): 
  print_unbuffered('')
  print_string_as_info('********* testing dict_to_defaultdict function ****************')
  print_unbuffered('')

  print('This is what a new SingletonCounts object looks like: ')
  print(Singleton_Counts()) 
  print_unbuffered('')

  singleton_counts = {4: [0, 1, 2, 0, 0]}
  print('singleton_counts in dict format:')
  print(singleton_counts) 
  print_unbuffered('')

  singleton_counts = dict_to_defaultdict(singleton_counts)
  print('after conversion from dict to defaultdict:') 
  print(singleton_counts) 
  print_unbuffered('')

  non_existent_key = 10
  print(f"value for key ({non_existent_key}) that doesn't exist:") 
  print(singleton_counts[10]) 
  print_unbuffered('')

  SNV_count = 1
  conditioned_singleton_counts = np.array([0, 1])
  # THIS ALSO WORKS: conditioned_singleton_counts = [0, 1]
  print('SNV and singleton count data to add:')
  print(SNV_count, conditioned_singleton_counts)
  print_unbuffered('')

  singleton_counts[SNV_count] += conditioned_singleton_counts  
  print('After adding SNV and singleton count data:')
  print(singleton_counts)
  print_unbuffered('')

  singleton_counts[SNV_count] += conditioned_singleton_counts  
  print('After adding SNV and singleton count data a second time:')
  print(singleton_counts)
  print_unbuffered('')

  SNV_count = 1
  conditioned_singleton_counts = np.array([2, 3])
  print('SNV and singleton count data to add:')
  print(SNV_count, conditioned_singleton_counts)
  print_unbuffered('')

  singleton_counts[SNV_count] += conditioned_singleton_counts  
  print('After adding SNV and singleton count data:')
  print(singleton_counts)
  print_unbuffered('')

  singleton_count = 1 
  singleton_counts[SNV_count][singleton_count] += 1
  print(f'After adding count data ({SNV_count}, {singleton_count}) from a single window...')
  print(singleton_counts)
  print_unbuffered('')

  print('singleton_counts items: ')
  print(singleton_counts.items())
  print_unbuffered('')

  print('types of values of singleton_counts: ')
  print([type(v) for v in singleton_counts.values()])
  print_unbuffered('')


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
  test_dict_to_defaultdict()
  test_data_structure()
  test_combine_singleton_counts()


