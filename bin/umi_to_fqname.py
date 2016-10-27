 #!/usr/bin/env python3
import argparse
import os, sys
import gzip
from itertools import islice
from collections import Counter

"""
Tools to extract UMIs from reads and append to another fastq file
Heavily borrowed from umitools https://github.com/brwnj/umitools

"""
IUPAC = {
    "A": "A",
    "T": "T",
    "C": "C",
    "G": "G",
    "R": "GA",
    "Y": "TC",
    "M": "AC",
    "K": "GT",
    "S": "GC",
    "W": "AT",
    "H": "ACT",
    "B": "GTC",
    "V": "GCA",
    "D": "GAT",
    "N": "GATC"
}

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
        
def get_umis(filehandle, umi):
    """extract read_ids:read_ids_UMIs line by line from fastq"""
  
    id_dict = {}
    umi_length = len(umi)
    umi_stats = Counter()
    for read in fq_reader(filehandle):
        read_name, umi_to_test = clip_umi(read, umi_length)
        if valid_umi(umi, umi_to_test):
            id_dict[read_name] = "{}:UMI_{}".format(read_name, umi_to_test)
        else:
            umi_stats.update([umi_to_test])
            continue
    
    print("Invalid UMI Total:", sum(umi_stats.values()), file=sys.stderr)
    print("Unique UMIs Removed:", len(list(umi_stats)), file=sys.stderr)
    print("Top", 50, "Invalid UMIs:", file=sys.stderr)
    for umi, val in umi_stats.most_common(50):
        print(umi, val, sep="\t", file=sys.stderr)
        
    return id_dict

def clip_umi(record, n):
    """Removes UMI sequence from read, and stores in a dictionary with
    read.
    Args:
        record (Fastq): `Fastq` record
        n (int): Length of the UMI, assumed to be on the 3' end       
    Returns:
        Fastq else str: The record or the failed UMI sequence
    >>> fq = Fastq(["@cluster_455 2","GGGGGAGCCACGAGGTGTGTTTTATTTTCATTATTC","+","C===>=B=@:<;4A;8=9?6EEC0?DDA72B@3EB4"])
    >>> r, umi = clip_umi(fq, 6,)
    >>> r
    Fastq(cluster_455:UMI_GGGGGA 2)
    """
    umi = record.seq[-n:]
    record.name = record.name.split(" ", 1)[0]
    
    return record.name, umi

def valid_umi(iupac, umi):
    """Parse UMI sequence to validate against IUPAC sequence.
    Args:
        iupac (str): IUPAC sequence
        umi (str): observed sequence
    Returns:
        bool
    >>> valid_umi("NNNV", "ACGT")
    False
    >>> valid_umi("NNNV", "ACGG")
    True
    """
    for code, base in zip(iupac, umi):
        try:
            if base not in IUPAC[code]:
                return False
        except KeyError:
            return False
    return True
 
def process_fastq(filehandle, umi_dict):  
    """ Parse a fastq and modify read name if read found in read:umi
    dictionary
    """
    
    for read2 in fq_reader(filehandle):
        read2_id = read2.name.split(" ", 1)[0]
        new_read2 = umi_dict.pop(read2_id, None)
        if new_read2:
            read2.name = new_read2
            print(read2)
        else:
            continue

def main():
    
    parser = argparse.ArgumentParser(description="""Add UMI information to a fastq,
                                       note that only the first field is used for
                                       pattern matching, other info is ignored (i.e. index info)""")
  
    parser.add_argument('-f1',
                          '--fastq_umi',
                          help ='fastq with UMI',
                       required = True)
    parser.add_argument('-f2',
                          '--fastq_mate',
                          help ='paired fastq to add UMI to read name',
                          required = True)
    parser.add_argument('-u',
                          '--umi_seq',
                          help ="""IUPAC representation of UMI, assuming that it
                          is at the 3p end, e.g. NNNN is a 4nt UMI at the 3p
                          end of read""",
                           required = True)
    args=parser.parse_args()
     
    f1 = args.fastq_umi
    f2 = args.fastq_mate
    umi = args.umi_seq
     
    gzopen = lambda f: gzip.open(f, 'rt') if f.endswith(".gz") else open(f)    

    f1 = gzopen(f1)
    f2 = gzopen(f2)
    umi_ids = get_umis(f1, umi)
    process_fastq(f2, umi_ids)

if __name__ == '__main__': main()



