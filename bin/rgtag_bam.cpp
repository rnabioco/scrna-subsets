#include <zlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "htslib/cram/sam_header.h"

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

BamReader::BamReader(const std::string& bampath, int cache_size) {
  const char* cbampath = bampath.c_str();
  in = sam_open(cbampath, "rb");
  if (in == NULL) {
    throw std::runtime_error("Fail to open BAM file") ;
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

std::unordered_set<std::string> get_rg_tags(std::string filename){
  
  BamReader bfile(filename) ;
  bam_hdr_t *h = sam_hdr_read(bfile.in) ; 
  
  std::unordered_set<std::string> rg_tgs ; 
  
  bam1_t *aln = bam_init1() ; 
  while (bam_read1(bfile.bz, aln) > 0) {
    uint8_t *tag_int = bam_aux_get(aln, "CN") ;
    char *tag_val = bam_aux2Z(tag_int ) ; 
    
    // convert to string to not deal with memory issues
    int len = strlen(tag_val) ;
    std::string s(tag_val, tag_val + len) ;
    rg_tgs.insert(s) ;
    
  }
  bam_destroy1(aln) ;
  
  return rg_tgs ; 
}

void write_new_bam(std::string oldbam, 
                   std::string newbam,
                   std::unordered_set<std::string>& rgtags){
  
  BamReader bfile(oldbam) ;
  bam_hdr_t *h = sam_hdr_read(bfile.in) ; // get header

  const char *fn_out = newbam.c_str() ;
  samFile *fp_out = sam_open(fn_out, "wb") ; //initialize output bam
  
  // add rg tags
  SAM_hdr *sh = sam_hdr_parse_(h->text, h->l_text) ;// get header for editing
  if (sh){
    for (auto i:rgtags){
      int l = sam_hdr_add(sh, "RG", "ID", i.c_str(), "SM", "UNKNOWN", NULL) ;
      if (l < 0) {
        throw std::runtime_error("error adding rg tag ") ;
      } 
    }
  }
  
  // write updated header
  if (int l = sam_hdr_rebuild(sh) < 0 ) {
    throw std::runtime_error("error updating header ") ;
  }
  // do a bunch of C stuff (taken from samtools/bam_reheader.c)
  free(h->text);
  h->text = strdup(sam_hdr_str(sh)); // update original header
  h->l_text = sam_hdr_length(sh); // update header lengths
  
  sam_hdr_free(sh) ;
  if (int l = sam_hdr_write(fp_out, h) < 0) {
    throw std::runtime_error("error writing header ") ;
  } 
  
  bam1_t *aln = bam_init1() ;
  while (bam_read1(bfile.bz, aln) > 0){
    uint8_t *tag_int = bam_aux_get(aln, "CN") ;
    char *tag_val = bam_aux2Z(tag_int) ; 
    
    // append prior to delete as *tag_val will change after deleting
    // tag_int
    bam_aux_append(aln, "RG", 'Z', strlen(tag_val) + 1, (uint8_t *)(tag_val)) ;
    bam_aux_del(aln, tag_int) ;
    
    int w = sam_write1(fp_out, h, aln) ;
    if (w < 0) {
      throw std::runtime_error("error writing alignment") ;
    } 
  }
  
  sam_close(fp_out) ;
}

int usage(){
  std::cerr << "./rgtag_bam <in.bam> <out.bam>" << std::endl ;
  return 1 ;
}

int main(int argc, char *argv[]){
  if (argc == 1) {
    return usage() ;
  } else if  (argc != 3) {
    std::cerr << "in bam and out bam arguments required" << std::endl ;
    return usage() ;
  }
  
  std::string filename(argv[1]) ;
  std::string outname(argv[2]) ;
  
  try {
    auto rg_tgs = get_rg_tags(filename) ;
    write_new_bam(filename, outname, rg_tgs) ;
  } catch (const std::runtime_error& e) {
    std::cerr << e.what();
  }

}
