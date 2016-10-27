#!/usr/bin/env python3
import argparse
import os, sys
import gzip
import pysam

def get_ids(ids):
  """extract read_ids:read_ids_UMIs line by line from gzipped text file"""
  
  ids = gzip.open(ids, 'rt')
  id_dict = {}
  for line in ids:
    read = line.split(" ")[0].strip("\n")[1:]
    read_id = read.split("UMI_")[0][:-1]
    id_dict[read_id] = read

  return id_dict
  ids.close()

def process_bam(ibam, obam, ids):
  """ iterate through bam file and rename fastq_id to 
     add UMI information
  """
  id_dict = get_ids(ids) 

  aln_file = pysam.AlignmentFile(ibam, "rb")
  out_file = pysam.AlignmentFile(obam, "wb", template = aln_file)

  for aln in aln_file.fetch(until_eof=True):
    read = aln.query_name
    read_umi = id_dict.pop(read, None)
    if read_umi:
      aln.query_name = read_umi
      out_file.write(aln)
    else:
      continue

  aln_file.close()
  out_file.close()

def main():
    
  parser = argparse.ArgumentParser(description="""Add UMI information to a fastq idline in bamfile,
                                   note that only the first field is used for
                                   pattern matching, other info is ignored (i.e. index info)""")
  
  parser.add_argument('-i',
                      '--ids',
                      help ='fastq ids with UMI, one per line',
                      required = True)
  
  parser.add_argument('-b',
                      '--bam',
                      help ='bam file',
                      required = True)
  parser.add_argument('-o',
                      '--output_file',
                      help = 'output file',
                      required = True)
  args=parser.parse_args()
  
  ids = args.ids
  bam = args.bam
  output_file = args.output_file
  
  process_bam(bam, output_file, ids)
  
    
if __name__ == '__main__': main()



