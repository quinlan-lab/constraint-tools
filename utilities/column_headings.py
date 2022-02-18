import gzip

def fetch_column_headings_indices(filename): 
  with gzip.open(filename, 'r') as fh: 
    column_headings = [b.decode('UTF-8') for b in fh.readline().strip().split()]
    heading_to_index = {
      column_heading: i
      for i, column_heading in enumerate(column_headings)
    }
    return column_headings, heading_to_index


