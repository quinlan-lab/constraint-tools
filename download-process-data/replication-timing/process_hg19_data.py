import sys 
import numpy as np 

def process_hg19_data(): 
  header = sys.stdin.readline()   
  header = f'#chromosome\tstart\tend\treplication timing per individual'
  header = header.strip()
  print(header)

  for line in sys.stdin: 
    fields = line.split('\t')

    region = '\t'.join(fields[:3])
    region = 'chr' + region

    replication_timing_all_individuals = []
    for field in fields[4:]: 
      try: 
        replication_timing_all_individuals.append(float(field))
      except ValueError: 
        continue     
    replication_timing_per_individual = np.mean(replication_timing_all_individuals)

    output = f'{region}\t{replication_timing_per_individual}'
    output = output.strip()
    print(output)

if __name__ == '__main__': 
  process_hg19_data()