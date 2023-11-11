import sys

def create_bed_file_of_large_windows(csv_filename):
  with open(csv_filename) as f: 
    header = f.readline() 
    for line in f: 
      line = '\t'.join(line.strip().split(','))
      print(line)

if __name__ == '__main__':
  create_bed_file_of_large_windows(csv_filename=sys.argv[1])