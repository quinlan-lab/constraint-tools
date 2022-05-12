# Use biomart for "Complex cross-database queries"
# Intro to biomart: 
# https://uswest.ensembl.org/info/data/biomart/how_to_use_biomart.html
# Martservices usage: 
# http://www.ensembl.org/biomart/martservice
# Obtaining the BioMart xml (used below) from the BioMart website: 
# https://uswest.ensembl.org/info/data/biomart/biomart_restful.html#biomartxml

import requests
import sys

try: 
  from colorize import print_string_as_info
except: 
  print_string_as_info = lambda text: print(text, file=sys.stderr) 

with open(sys.argv[1]) as f:
  xml = f.read()

with (
  # https://docs.python-requests.org/en/latest/user/advanced/#body-content-workflow
  requests.post('http://www.ensembl.org/biomart/martservice', data={'query': xml}, stream=True) as r
):
  # https://docs.python-requests.org/en/latest/api/#requests.Response.iter_lines
  # https://docs.python-requests.org/en/latest/user/advanced/#streaming-requests
  for line in r.iter_lines(decode_unicode=True):
    if line == '[success]':
      print_string_as_info('BioMart download is complete.')  
      break  
    print(line, file=sys.stdout)  

