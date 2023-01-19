import sys 
import numpy as np 

def process_hg19_data(): 
  for line in sys.stdin: 
    fields = line.split('\t')

    chromosome, position = fields[:2] 
    position = int(position)
    start = position - 1 
    end = start + 1 
    region = '\t'.join([chromosome, str(start), str(end)])

    ID, REF, ALT = fields[2:]

    if ID.startswith('1000G'): continue 

    ALT = ALT.strip()
    if (len(REF) > 1) or (len(ALT) > 1): 
      print(ID, REF, ALT)
      raise ValueError

    output = f'{region}\t{ID}\t{REF}\t{ALT}'
    print(output)

if __name__ == '__main__': 
  process_hg19_data()
