import sys 

with (
  open(sys.argv[1]) as data_hg38_file
): 
  header = 'chromosome\tstart\tend\tclinvar_id\tREF\tALT'
  print(header)

  for line in data_hg38_file: 
    print(line, end='')



  