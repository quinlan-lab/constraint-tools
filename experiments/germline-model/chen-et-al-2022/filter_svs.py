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
  with open('/scratch/ucgd/lustre-work/quinlan/u0055382/SV_constraint/SV_data/CCDG_DEL_singleton_ids.txt') as fh: 
    CCDG_singletons = set(line.strip() for line in fh)

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
    if source != sys.argv[1]: continue 

    if source == 'CCDG':
      alt_allele_count = 1 if sv_id in CCDG_singletons else '>1'
      number_individuals_who_were_genotyped = '.'
    else:
      alt_allele_count = 1*int(number_het_individuals) + 2*int(number_homalt_individuals)
      number_individuals_who_were_genotyped = (
        int(number_homref_individuals) + 
        int(number_het_individuals) + 
        int(number_homalt_individuals)
      )
      if number_individuals_who_were_genotyped < int(sys.argv[2]): continue 

    line = '\t'.join([
      chromosome, 
      start, 
      end, 
      sv_length, 
      sv_id, 
      str(alt_allele_count), 
      str(number_individuals_who_were_genotyped)
    ])
    print(line)

if __name__ == '__main__': 
  print_head() 
  try:
    print_tail()
  except BrokenPipeError: 
    pass