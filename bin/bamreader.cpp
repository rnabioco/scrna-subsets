#include "subsampling.h"

BamReader::BamReader(const std::string& bampath, 
                     bool check_idx,
                     int cache_size) {
  const char* cbampath = bampath.c_str();
  in = sam_open(cbampath, "rb");
  if (in == NULL) {
    throw std::runtime_error("Fail to open BAM file") ;
  }
  
  bz = in->fp.bgzf ; // bgzf file pointer
  idx = bam_index_load(cbampath); // load BAM index
  if (idx == 0 && check_idx) {
    std::cerr << "BAM indexing file is not available for file " << bampath << std::endl;
  }
  
  if (cache_size > 0){
    bgzf_set_cache_size(bz, cache_size);
  }
}