#!/usr/bin/env python3
import argparse
import os, sys
import gzip

def extract_ids(fq_in, ids_out):
  """ extract out read_ids to map for demultiplexing R2
  """  
 
  fq = gzip.open(fq_in, 'rt')

  with gzip.open(ids_out, 'at') as in_file:
    while True:
  #parse fastq
      idline = fq.readline()
      if not idline:
        break
      seq    = fq.readline()
      spacer = fq.readline()
      quals  = fq.readline()
      read_id = idline.split()[0].strip("@")
      file_id = str(fq_in)
      data = '{}\n'.format( read_id)
      in_file.write(data) 

def main():

  parser = argparse.ArgumentParser(description="""extract read_ids from fastq.gz""")

  parser.add_argument('-f',
                      '--fastq',
                      help ='fastq gzipped',
                      required = True)

  parser.add_argument('-o',
                      '--output',
                      help = 'output txt.gz',
                      required = True)
  args=parser.parse_args()

  fq = args.fastq
  output = args.output

  extract_ids(fq, output)


if __name__ == '__main__': main()
