import pyranges as pr

# https://biocore-ntnu.github.io/pyranges/loadingcreating-pyranges.html

def get_all_exons(): 
  return pr.read_bed('/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/genes/grch38/exons.sorted.bed.gz')

def get_constitutive_exons(): 
  return pr.read_bed('/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/genes/grch38/constitutive-exons.sorted.bed.gz')

def get_canonical_exons(): 
  return pr.read_bed('/scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/genes/grch38/canonical-exons.sorted.bed.gz')
  
if __name__ == '__main__': 
  print(get_constitutive_exons()['chr20', 63468750:63473750])