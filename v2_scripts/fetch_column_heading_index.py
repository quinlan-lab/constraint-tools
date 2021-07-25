import sys 
from column_headings import fetch_column_heading_index

if __name__ == '__main__': 
  print(fetch_column_heading_index(filename=sys.argv[1], column_heading=sys.argv[2], unit_offset=True))
