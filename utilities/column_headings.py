import gzip

def fetch_column_heading_indices(filename): 
  with gzip.open(filename, 'r') as fh: 
    column_headings = [b.decode('UTF-8').lower() for b in fh.readline().strip().split()]
  return {
    column_heading: i
    for i, column_heading in enumerate(column_headings)
  }

def fetch_column_heading_index(filename, column_heading, unit_offset=False): 
  column_heading_indices = fetch_column_heading_indices(filename)
  column_heading_index = column_heading_indices[column_heading.lower()]
  if unit_offset: 
    return column_heading_index + 1
  else:
    return column_heading_index


