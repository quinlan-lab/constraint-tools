import sys 

from pack_unpack import unpack 

target_label = sys.argv[1] # "positive" or "negative"

def convert_to_bed(): 
  for line in sys.stdin: 
    fields = line.split('|')
    region = fields[1].strip()
    element_label = fields[3].strip()
    if element_label != target_label: continue 
    chromosome, start, end = unpack(region)
    print(f'{chromosome}\t{start}\t{end}')

if __name__ == '__main__': 
  convert_to_bed()