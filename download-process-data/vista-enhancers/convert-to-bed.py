import sys 

from pack_unpack import unpack 

def convert_to_bed(): 
  for line in sys.stdin: 
    region = line.split('|')[1].strip()
    chromosome, start, end = unpack(region)
    print(f'{chromosome}\t{start}\t{end}')

if __name__ == '__main__': 
  convert_to_bed()