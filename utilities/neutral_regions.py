import pyranges as pr

from pack_unpack import unpack 

def get_all_neutral_regions(model): 
  # https://biocore-ntnu.github.io/pyranges/loadingcreating-pyranges.html
  return pr.read_bed(model['neutralRegions'])

def get_neutral_regions(region, model): 
  neutral_regions = get_all_neutral_regions(model)
  
  chromosome, start, end = unpack(region)
  # https://biocore-ntnu.github.io/pyranges/manipulating-the-data-in-pyranges.html
  return neutral_regions[chromosome, start:end].df.to_dict(orient='records')

