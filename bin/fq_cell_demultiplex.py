#!/usr/bin/env python

import sys
import os
import argparse
from collections import defaultdict

def seq_differences(str_a, str_b):
  """ enumerate the number of mismatches between two strings
  """
  return(sum(1 for x,y in zip(str_a, str_b) if x != y))

def demultiplex_known(fq1, fq2, directory, input_barcodes, mismatch_allowance):
  """ identify cell barcodes from read 1 and generate dictionary to
    filter fastq with.
    fastq parsing idea https://www.biostars.org/p/10353/
  """
  fq1 = open(fq1, 'r')
  barcode_ids = defaultdict(list)
  
  known_barcodes = set( x.strip() for x in input_barcodes )
  
  while True:
    idline = fq1.readline()
    if not idline:
      break
    seq    = fq1.readline()
    spacer = fq1.readline()
    quals  = fq1.readline()
    barcode = seq.strip("\n")[:11]
    read_id = idline.strip()[:-2]
    if barcode in known_barcodes:
      barcode_ids[barcode].append(read_id)
      
    elif mismatch_allowance == 0:
        continue
    
    else:
      scores = []
      for test_barcode in known_barcodes:
        scores.append([barcode, test_barcode, seq_differences(barcode, test_barcode)])
      pass_filter = []
      for alignment in scores:
        if alignment[2] <= int(mismatch_allowance):
          pass_filter.append(alignment)
      if len(pass_filter) != 1:
        continue
      else:
        barcode_ids[pass_filter[0][1]].append(read_id)
  print("hello")
  if len(barcode_ids.values()) > 0:
    split_fastq(fq2, barcode_ids, directory)
  else:
    sys.exit("no barcodes found in %s" % fq1.name())
    
def demultiplex_unknown(fq1, fq2, directory):
  """generate dictionary of barcodes and reads Ids from read 1
     then parse read 2 to extract entries into distinct files
  """
  fq1 = open(fq1, 'r')
  barcode_ids = {}
  
  while True:
    idline = fq1.readline()
    if not idline:
      break
    seq    = fq1.readline()
    spacer = fq1.readline()
    quals  = fq1.readline()
    barcode = seq.strip("\n")[:11]
    read_id = idline.strip()[:-2]
    barcode_ids[read_id] = barcode
    
  if len(barcode_ids.values()) > 0:
    split_fastq(fq2, barcode_ids, directory)
  else:
    sys.exit("no barcodes found in %s" % fq1.name())
    

def split_fastq(fq, filter_dict, directory):
  """ split fastq records into new files with barcode as name
     filter_dict key:values are read_id:barcode
     
  """
  directory = directory.strip("/")
  os.system('mkdir -p ./' + directory)
  os.system('rm -f ./' + directory + "/*")
  
  fq = open(fq, 'r') 
  while True:
  #parse fastq
    idline = fq.readline()
    if not idline:
      break
    seq    = fq.readline()
    spacer = fq.readline()
    quals  = fq.readline()
    read_id = idline.strip()[:-2]
    
    #if filter_dict[read_id]:
    barcode = filter_dict[read_id]
    filename = directory + "/" + "Cell_" + barcode + ".fastq"
    filename = open(filename, 'a')
    filename.write('%s%s%s%s' % ( idline, seq, spacer, quals ))
  filename.close()

  
   #report barcode counts to stderr
  
  counts = {}
  for read, barcode in filter_dict.iteritems():
    if barcode in counts:
      counts[barcode] += 1
    else:
      counts[barcode] = 1
            
  sys.stderr.write('barcode_sequence\tcount\n')
  for w in sorted(counts, key=counts.get, reverse=True):
    sys.stderr.write('%s\t%s\n' % (w, counts[w]))
    
def main():
    
  parser = argparse.ArgumentParser(description="""Demultiplex paired-end Wafergen-Style fastq reads""")
  
  parser.add_argument('-f1',
                      '--left_mate_fastq',
                      help ='Read 1 containing the cell barcode and the UMI',
                      required = True)
  
  parser.add_argument('-f2',
                      '--right_mate_fastq',
                      help = 'Read 2 containing the cDNA read',
                      required = True)
  parser.add_argument('-a',
                      '--barcode_seqs',
                      help = 'tsv file with barcodes listed as ID tab Barcode, one per line',
                      required = False)
  parser.add_argument('-d',
                      '--directory',
                      help = 'output directory',
                      required = True)
  parser.add_argument('-m',
                      '--mismatches',
                      help = 'int number of mismatches allowed between read barcode and known barcode, defaults to 1',
                      required = False)
  
  args=parser.parse_args()
  
  fq_1 = args.left_mate_fastq
  fq_2 = args.right_mate_fastq
  directory = args.directory
  mismatches = args.mismatches
  
  if not mismatches:
    mismatches = 1
    
  if args.barcode_seqs:
    input_barcodes = open(args.barcode_seqs, 'r')
    demultiplex_known(fq_1, fq_2, directory, input_barcodes, mismatches)
  else:
    demultiplex_unknown(fq_1, fq_2, directory)
    
    
if __name__ == '__main__': main()
  