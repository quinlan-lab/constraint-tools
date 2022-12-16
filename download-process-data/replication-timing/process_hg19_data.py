import sys 

def process_hg19_data(): 
  header = sys.stdin.readline()   
  replication_timing_header = '\t'.join(header.split('\t')[4:])
  header = f'#chromosome\tstart\tend\t{replication_timing_header}'
  header = header.strip()
  print(header)

  for line in sys.stdin: 
    fields = line.split('\t')
    region = '\t'.join(fields[:3])
    region = 'chr' + region
    replication_timing_data = '\t'.join(fields[4:])
    output = f'{region}\t{replication_timing_data}'
    output = output.strip()
    print(output)

if __name__ == '__main__': 
  process_hg19_data()