import sys 

def print_head(): 
  header = sys.stdin.readline() 
  header = '\t'.join([
    'chromosome', 
    'start', 
    'end', 
    'sv_length', 
    'sv_id', 
    'alt_allele_count', 
    'number_individuals_who_were_genotyped'
    ])
  print(header)

def print_tail(): 
  for line in sys.stdin: 
    fields = line.strip('\n').split('\t')
    (
      chromosome, 
      start, 
      end, 
      sv_length, 
      sv_type, 
      source, 
      sv_id, 
      allele_frequency, 
      number_homref_individuals, 
      number_het_individuals, 
      number_homalt_individuals
    ) = fields[:11]
    chromosome = f'chr{chromosome}'
    if int(sv_length) < 100: continue 
    if sv_type != 'DEL': continue 
    if source != 'gnomAD': continue 
    alt_allele_count = number_het_individuals  
    number_individuals_who_were_genotyped = (
      int(number_homref_individuals) + 
      int(number_het_individuals) + 
      int(number_homalt_individuals)
    )
    if number_individuals_who_were_genotyped < 10000: continue 
    line = '\t'.join([
      chromosome, 
      start, 
      end, 
      sv_length, 
      sv_id, 
      alt_allele_count, 
      str(number_individuals_who_were_genotyped)
    ])
    print(line)

if __name__ == '__main__': 
  print_head() 
  try:
    print_tail()
  except BrokenPipeError: 
    pass