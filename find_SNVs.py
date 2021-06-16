# pysam API: 
# https://github.com/pysam-developers/pysam/blob/b82cbcae22c088e64fdb58f8acaf1e9773c7b088/pysam/libctabix.pyx
import pysam

mutations = pysam.TabixFile("data/icgc/mutations.sorted.maf.gz")

print(mutations.contigs)

for row in mutations.fetch("1", 2000000, 2010000, parser=pysam.asTuple()):
  print('*************************')
  print(type(row))
  print(row[0])
  print(row)
