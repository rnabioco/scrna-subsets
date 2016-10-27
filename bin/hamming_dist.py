#!/usr/bin/env python
import gzip
import sys
import os
import argparse
from collections import defaultdict, Counter
from umi_to_fqname import fq_reader, Fastq

def seq_differences(str_a, str_b):
  """ enumerate the number of mismatches between two strings
  """
  return(sum(1 for x,y in zip(str_a, str_b) if x != y))

def extract_bc(in_file):
  bcs = []
  for line in open(in_file):
    bcs.append(line.split("\t")[1])
  return(bcs)
  in_file.close()

def main():
    
  parser = argparse.ArgumentParser(description="""calculate minimum hamming distance between all provided sequences """)
  
  parser.add_argument('-i',
                      '--input',
                      help ='barcodes listed in second column',
                      required = True)
  parser.add_argument('-f',
                      '--fastq',
                      help = "fastq file to extract barcodes from",
                      required = True)

  gzopen = lambda f: gzip.open(f, 'rt') if f.endswith(".gz") else open(f)
  
  args=parser.parse_args()
  
  bcs = args.input
  fq = args.fastq  
  fq = gzopen(fq)

  bcs = extract_bc(bcs)
  cnt = Counter()
  for read in fq_reader(fq):
    bc_hamming_each = []
    read.seq = read.seq[:11]
    for barcode in bcs:
      bc_diff = seq_differences(read.seq, barcode)
      bc_hamming_each.append(bc_diff)
    bc_min = min(bc_hamming_each)
    cnt.update([bc_min])

  print("hamming_dist\tcounts")
  for hamming_dist, counts in cnt.most_common():
    print("{}\t{}".format(hamming_dist, counts))
    
if __name__ == '__main__': main()
  
