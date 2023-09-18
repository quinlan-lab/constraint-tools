import sys 

from pack_unpack import unpack 

def remove_prefix(text, prefix):
  if text.startswith(prefix):
    return text[len(prefix):]
  return text

def convert_to_bed(): 
  target_class = sys.argv[1] # "positive" or "negative"

  count = 0
  element_id = None
  for line in sys.stdin: 
    line = remove_prefix(line, '<pre>') # curl pulls down the <pre> HTML tag, so remove it

    if not line.startswith('>Human'): continue

    fields = line.split('|')
    region = fields[1].strip()
    element_id = fields[2].strip()
    element_class = fields[3].strip()

    # The element ids are not consecutive, implying that 
    # the number of vista elements in the webpage is smaller than the id of the last vista element.
    count += 1
    element_id = int(remove_prefix(element_id, 'element').strip())
    print(count, element_id, file=sys.stderr)

    if element_class != target_class: continue 
    chromosome, start, end = unpack(region)
    print(f'{chromosome}\t{start}\t{end}')

  # I count 1947 human vista elements, but the id of the last human vista element in the webpage is 2668.
  print(f'Number of human vista enhancers: {count}', file=sys.stderr)
  print(f'ID of last human vista element: {element_id}', file=sys.stderr)

if __name__ == '__main__': 
  convert_to_bed()