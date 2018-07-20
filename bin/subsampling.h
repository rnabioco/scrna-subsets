#ifndef subsampling__subsampling_H 
#define subsampling__subsampling_H

#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <vector>
#include <string>
#include <map>
#include <unordered_set>
#include <set>
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"

// class for handling bam file opening and closing
class BamReader {
public:
  samFile* in;
  hts_idx_t* idx;
  BGZF* bz;
  BamReader(const std::string& bampath, 
            bool check_idx = true,
            int cache_size=10*BGZF_MAX_BLOCK_SIZE) ;
  
  ~BamReader(){
    hts_idx_destroy(idx);
    sam_close(in);
  }
};

// Class for handling umis
class Molecule {
public:
  std::string cn, bo, xt, umikey ;
  int64_t ug ; 
  bool has_tags ; 
  Molecule(bam1_t*) ;
};

// Class for handling umis with positions
class PosMolecule {
public:
  int32_t pos ;
  std::string chrom, cn, bo, xt, umikey ;
  bool has_tags ; 
  PosMolecule(bam1_t*, bam_hdr_t*);
};

#endif
