#!/usr/bin/env python
import argparse
import os, sys
import gzip
from itertools import chain, islice

def chunks(iterable, n):
   "chunks(ABCDE,2) => AB CD E"
   iterable = iter(iterable)
   while True:
       yield chain([next(iterable)], islice(iterable, n-1))


def split_fastq(fq_in, line_count, suffix, prefix):
  line_count = int(line_count) 
  with gzip.open(fq_in, 'rb') as bigfile:
      for idx, lines in enumerate(chunks(bigfile, line_count)):
          file_split = '{}_{}{}'.format(prefix, idx, suffix)
          with gzip.open(file_split, 'wb') as f:
              f.writelines(lines)


def main():

  parser = argparse.ArgumentParser(description="""extract read_ids and filenames from fastq.gz""")

  parser.add_argument('-f',
                      '--fastq',
                      help ='fastq gzipped',
                      required = True)

  parser.add_argument('-l',
                      '--line_count',
                      help = 'number of lines to chunk to each new file',
                      required = True)
  parser.add_argument('-s',
                      '--suffix',
                      help = 'suffix for output filename', 
                      required = True)
  parser.add_argument('-p',
                      '--prefix',
                      help = 'prefix for output filename',
                      required = True)
  args=parser.parse_args()

  fq = args.fastq
  line_count = args.line_count
  suffix = args.suffix
  prefix = args.prefix

  split_fastq(fq, line_count, suffix, prefix)


if __name__ == '__main__': main()
