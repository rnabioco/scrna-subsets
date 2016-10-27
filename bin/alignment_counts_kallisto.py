#!/usr/env/bin python3

import argparse
import os, sys
import gzip
import pysam 
import pdb

def kallisto_map_stats(pbam_in, stat_out):
  ''' calculate total reads, single aligned, multiple aligned, unaligned from bam '''

  track = pysam.AlignmentFile(pbam_in, "rb")
  
  total_reads = set()
  unique_aln = 0
  multiple_aln = set()
  unmapped_reads = 0
  
  for aln in track.fetch(until_eof=True):
    flag = aln.flag
    read = aln.query_name
    try:
      tag = aln.get_tag("NH")
    except KeyError:
      tag = 0
    
    if flag == 4:
      unmapped_reads += 1
    elif flag == 0 and tag == 1:
      unique_aln += 1 
    else:
      multiple_aln.add(read)

    total_reads.add(read)
  
  multiple_aln = len(multiple_aln)
  total_reads = len(total_reads)
  try:
    percent_aln = 100 * (unique_aln + multiple_aln) / total_reads
  except ZeroDivisionError:
    percent_aln = 0
  percent_unaln = 100 - percent_aln

  out = open(stat_out, 'wt')
  file_id = str(pbam_in)
  file_id = os.path.basename(file_id)

  stats = [total_reads, unique_aln, multiple_aln, unmapped_reads, percent_aln, percent_unaln]
  stat_name = ["total reads", "unique alignments", "multiple alignments", "unmapped reads", "percent alignment", "percent unmapped"]
  
  for idx, stat in enumerate(stats):
    out.write("{}\t{}\t{}\n".format(file_id, stat_name[idx], stat))


def main():

  parser = argparse.ArgumentParser(description="""calculate alignment stats from kallisto pseudobam""")
  
  parser.add_argument('-i', '--bam', help = 'input pseudobam file from kallisto', required = True)

  parser.add_argument('-o', '--output', help = 'output text file', required = True)
  
  args=parser.parse_args()

  pbam_in = args.bam
  stat_out = args.output

  kallisto_map_stats(pbam_in, stat_out)
  
if __name__ == '__main__': main ()

