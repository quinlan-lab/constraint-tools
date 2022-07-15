import pyranges as pr

from pack_unpack import unpack 

def get_all_trustworthy_noncoding_regions(trustworthy_noncoding_regions_filename): 
  # https://biocore-ntnu.github.io/pyranges/loadingcreating-pyranges.html
  return pr.read_bed(trustworthy_noncoding_regions_filename)

def get_trustworthy_noncoding_regions(region, trustworthy_noncoding_regions_filename): 
  trustworthy_noncoding_regions = get_all_trustworthy_noncoding_regions(trustworthy_noncoding_regions_filename)
  chromosome, start, end = unpack(region)
  # https://biocore-ntnu.github.io/pyranges/manipulating-the-data-in-pyranges.html
  return trustworthy_noncoding_regions[chromosome, start:end].df.to_dict(orient='records')

