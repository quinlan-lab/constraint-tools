def unpack(region): 
  chromosome, start_end = region.split(':')
  start, end = map(lambda s: int(s.replace(',', '')), start_end.split('-'))
  return chromosome, start, end

def pack(chromosome, start, end): 
  return '{}:{}-{}'.format(chromosome, start, end)

if __name__ == '__main__': 
  print(pack(1, 1, 100))