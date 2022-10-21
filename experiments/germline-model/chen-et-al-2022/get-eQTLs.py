#!/usr/bin/env python

import sys 

from colorize import print_string_as_info 

def get_eQTLs(): 
  with open(sys.argv[1], 'r') as fh:
    header = next(fh) # skip header

    for line in fh: 
      fields = line.split('\t')
      chromosome, position = fields[0].split(':')[:2]
      position = int(position)
      ws_posterior = fields[-1]
      print(f'chr{chromosome}\t{position-1}\t{position}\t{ws_posterior}')

if __name__ == '__main__': 
  get_eQTLs()