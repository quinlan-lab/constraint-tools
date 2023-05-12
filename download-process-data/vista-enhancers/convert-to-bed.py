import sys 

from pack_unpack import unpack 

def convert_to_bed(): 
  for line in sys.stdin: 
    fields = line.split('|')
    region = fields[1].strip()
    enhancer_class = fields[3].strip()
    if enhancer_class == 'negative': continue 
    chromosome, start, end = unpack(region)
    print(f'{chromosome}\t{start}\t{end}')

if __name__ == '__main__': 
  convert_to_bed()