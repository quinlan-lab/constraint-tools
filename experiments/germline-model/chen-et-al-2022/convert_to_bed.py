import sys 

for line in sys.stdin: 
  line = line.strip('\n')
  fields = line.split('\t')
  chromosome, position = fields[:2]
  other_fields = '\t'.join(fields[2:])
  position = int(position)
  start, end = position-1, position
  print(f'{chromosome}\t{start}\t{end}\t' + other_fields)