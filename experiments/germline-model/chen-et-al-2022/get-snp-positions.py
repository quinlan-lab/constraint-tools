#!/usr/bin/env python

import sys 

from colorize import print_string_as_info 

def get_snp_positions(): 
  with open(sys.argv[1], 'r') as fh:
    header = next(fh) # skip header

    count_empty = 0
    count_interaction = 0
    count_multiple = 0
    count_lines = 0 

    for line in fh: 
      count_lines += 1

      chromosome, position = line.split('\t')[11:13]

      if chromosome == '': 
        count_empty += 1
        continue 
      if 'x' in chromosome: 
        count_interaction +=1 
        continue         
      if ';' in chromosome:
        count_multiple += 1 
        continue 

      try: 
        position = int(position)
      except: # debugging
        for key, value in zip(header.split('\t'), line.split('\t')):
          print(key, value)
        sys.exit(1)

      print(f'chr{chromosome}\t{position-1}\t{position}')

    print_string_as_info(f'Fraction of GWAS-catalog records with empty chromosome and position fields: {count_empty/count_lines}')
    print_string_as_info(f'Fraction of GWAS-catalog records representing interacting SNPs: {count_interaction/count_lines}')
    print_string_as_info(f'Fraction of GWAS-catalog records reporting multiple SNPs: {count_multiple/count_lines}')

if __name__ == '__main__': 
  get_snp_positions()