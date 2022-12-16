import sys 

with (
  open(sys.argv[1]) as replication_timing_data_hg19_file,
  open(sys.argv[2]) as replication_timing_data_hg38_file
): 
  header = replication_timing_data_hg19_file.readline()   
  header = header.lstrip('#')
  print(header, end='')

  for line in replication_timing_data_hg38_file: 
    print(line, end='')



  