#!/usr/bin/env python3
import argparse
import os, sys
import gzip
import re
from itertools import islice
from collections import Counter
import pysam

"""
check read structure conforms to known wafergen style read 1 structure
([ATCG]{21}TTTTT
Class method borrowed from umitools https://github.com/brwnj/umitools

"""

class Fastq(object):
    """FASTQ record. Holds record of name, sequence, and quality scores.
    >>> fq = Fastq(["@illumina_naming_scheme", "ACTGACTG", "+", "KKKKKKKK"])
    >>> fq.name, fq.seq, fq.qual
    ('illumina_naming_scheme', 'ACTGACTG', 'KKKKKKKK')
    >>> fq
    Fastq(illumina_naming_scheme)
    >>> print(fq)
    @illumina_naming_scheme
    ACTGACTG
    +
    KKKKKKKK
    >>> fq = Fastq(["@fail", "ACTGACTG", "+", "KKK"]) # doctest: +ELLIPSIS
    Traceback (most recent call last):
     ...
    AssertionError: Seq and Qual vary in length
    """
    def __init__(self, args):
        self.name = args[0][1:]
        self.seq = args[1]
        self.qual = args[3]
        assert len(self.seq) == len(self.qual), "Seq and Qual vary in length"

    def __repr__(self):
        return "Fastq(%s)" % self.name

    def __str__(self):
        return "@%s\n%s\n+\n%s" % (self.name, self.seq, self.qual)

def fq_reader(filehandle):
    """Fastq iterator.
    Args:
        filehandle (file): open file handle
    Yields:
        Fastq
    """

    fqclean = (x.strip("\r\n") for x in filehandle if x.strip())
    while True:
        rd = [x for x in islice(fqclean, 4)]
        if not rd:
            raise StopIteration
        assert all(rd) and len(rd) == 4
        yield Fastq(rd) 
        
def examine_read_structure(filehandle, outfile):
  
    read_stats = Counter()
    read_stats["correct structure"] = 0
    read_stats["incorrect structure"] = 0

    read_regex = "[ATCGN]{21}TTTTT$"

    for read in fq_reader(filehandle):
      read_match = re.search(read_regex, read.seq)
      
      if read_match is not None:    
          read_stats.update(["correct structure"])
      else:
          read_stats.update(["incorrect structure"])
          
    for stat, val in read_stats.items():
        outfile.write("{}\t{}\n".format(stat, val))
        
def examine_bam_read_structure(bamfilehandle, outfile):
    '''extract umi from aligned bam file read name and compare to expected N{10}T{5} '''

    umi_stats = Counter()
    umi_stats["correct structure"] = 0
    umi_stats["incorrect structure"] = 0

    umi_regex = "[ATCGN]{10}TTTTT$"

    for read in bamfilehandle: 
     
      if read.is_unmapped:
          continue
      else :
          read_id = read.qname
          umi = read_id.split("_")[1]
          umi_match = re.search(umi_regex, umi)
          
          if umi_match is not None:
              print(umi)
              umi_stats.update(["correct structure"])
          else:
              umi_stats.update(["incorrect structure"])
          
    for stat, val in umi_stats.items():
        outfile.write("{}\t{}\n".format(stat, val))

def main():
    
  parser = argparse.ArgumentParser(description="""Check read structure of
          Wafergen-Style fastq reads\n
          Expect N{21}TTTTT """)
  
  parser.add_argument('-f',
                      '--fastq',
                      help ='Read 1 containing the cell barcode and the UMI',
                      required = False)
  parser.add_argument('-b',
                      '--bam',
                      help = "bam file input with UMI in read name",
                      required = False)
  parser.add_argument('-o',
                      '--out',
                      help = "output file",
                      required = True)
  
  args=parser.parse_args()
  
  fq = args.fastq
  bam = args.bam
  out = args.out
    
  gzopen = lambda f: gzip.open(f, 'rt') if f.endswith(".gz") else open(f)    

  if fq:
      fq = gzopen(fq)
      out = open(out, mode = 'w')
      examine_read_structure(fq, out)
      out.close()
  elif bam:
      try:
        samfile = pysam.AlignmentFile(bam, "rb")
      except ValueError:
        print("empty file", file = sys.stderr)  
        open(out, mode = 'a').close()
        sys.exit(0)   

      out = open(out, mode = 'w')
      examine_bam_read_structure(samfile, out)
      out.close()
  else:
      sys.exit("Either a bam or a fastq must be supplied")

if __name__ == '__main__': main()
  
