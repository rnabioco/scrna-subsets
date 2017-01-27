#include <zlib.h>
#include <stdio.h>
#include <regex>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <cstring>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/kseq.h>

// class for handling bam file opening and closing
class BamReader {
public:
  samFile* in;
  hts_idx_t* idx;
  BGZF* bz;
  BamReader(const std::string& bampath, int cache_size=10*BGZF_MAX_BLOCK_SIZE) ;
  
  ~BamReader(){
    hts_idx_destroy(idx);
    sam_close(in);
  }
};

// class for handling bam file opening and closing
BamReader::BamReader(const std::string& bampath, int cache_size) {
  const char* cbampath = bampath.c_str();
  in = sam_open(cbampath, "rb");
  if (in == NULL) {
    std::cerr << "Fail to open BAM file " << bampath << std::endl;
  }
  
  bz = in->fp.bgzf ; // bgzf file pointer
  idx = bam_index_load(cbampath); // load BAM index
  if (idx == 0) {
    std::cerr << "BAM indexing file is not available for file " << bampath << std::endl;
  }
  
  if (cache_size > 0){
    bgzf_set_cache_size(bz, cache_size);
  }
}


int main(int argc, char *argv[])
{
  if (argc == 1) {
    std::cerr << "Usage: %s <in.bam> \n" << argv[0] ;
    return 1;
  } 
  
  
  //read in bam
  BamReader bfile(argv[1]) ;  //handle opening bam
  const bam_hdr_t *header = sam_hdr_read(bfile.in) ; // get header
  bam1_t *aln = bam_init1() ; // initialize empty alignment container
  
  typedef std::map<std::string, int> barcodeDict ; 
  barcodeDict bcs ; 
  while (bam_read1(bfile.bz, aln) > 0) { // negative return values are errors
    std::string cn = bam_aux2Z(bam_aux_get(aln, "CN")) ;
    std::string cb = bam_aux2Z(bam_aux_get(aln, "CB")) ;
    std::string cbc = cn + "::" + cb ;
    bcs[cbc] += 1 ; 
  }
  
  for(auto const& it : bcs) {
    std::cout << it.first << "\t" << it.second << std::endl ; 
  }
  bam_destroy1(aln) ;

  return 0;
}
