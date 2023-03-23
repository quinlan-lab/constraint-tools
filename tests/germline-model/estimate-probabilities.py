# expected_singleton_counts = XXX

import sys 
model_filename = sys.argv[1]

import json 
with open(model_filename, 'r') as fh:
  d = json.load(fh) 
  observed_singleton_counts = d['singletonCounts']

from unittest import TestCase
# TestCase().assertDictEqual(observed_singleton_counts, expected_singleton_counts)

from colorize import print_string_as_info
# print_string_as_info('*** TESTS PASSED ***')



 
