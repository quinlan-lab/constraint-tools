expected_singleton_counts = {
  "6": [
    2,
    0,
    0,
    0,
    2,
    0,
    0
  ],
  "8": [
    0,
    2,
    0,
    0,
    2,
    0,
    2,
    0,
    0
  ],
  "5": [
    0,
    0,
    2,
    0,
    0,
    0
  ],
  "1": [
    2,
    2
  ],
  "2": [
    2,
    4,
    4
  ],
  "3": [
    2,
    4,
    2,
    0
  ],
  "7": [
    0,
    0,
    4,
    0,
    0,
    0,
    0,
    0
  ],
  "10": [
    0,
    0,
    0,
    0,
    0,
    2,
    2,
    0,
    0,
    0,
    0
  ],
  "4": [
    0,
    4,
    2,
    0,
    0
  ],
  "16": [
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    2,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0
  ],
  "9": [
    0,
    0,
    0,
    0,
    0,
    0,
    4,
    0,
    0,
    0
  ],
  "14": [
    0,
    0,
    0,
    0,
    0,
    0,
    2,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0
  ],
  "18": [
    0,
    0,
    0,
    0,
    0,
    2,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0
  ],
  "19": [
    0,
    0,
    0,
    0,
    0,
    0,
    2,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0
  ]
} 

import sys 
model_filename = sys.argv[1]

import json 
with open(model_filename, 'r') as fh:
  d = json.load(fh) 
  observed_singleton_counts = d['singletonCounts']

from unittest import TestCase
TestCase().assertDictEqual(observed_singleton_counts, expected_singleton_counts)

from colorize import print_string_as_info
print_string_as_info('*** TESTS PASSED ***')



 
