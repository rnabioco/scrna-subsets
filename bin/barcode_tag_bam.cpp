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

void splitName(const std::string &name, char delim, std::vector<std::string> &elems) {
    std::stringstream ss;
    ss.str(name);
    std::string item;
    while (std::getline(ss, item, delim)) {
      elems.push_back(item);
    }
  }

// leaves a trailing delim...
std::string join(const std::vector<std::string>& in_vector, const char* delim)
{
  std::stringstream res;
  copy(in_vector.begin(), in_vector.end(), std::ostream_iterator<std::string>(res, delim));
  return res.str();
}

std::map<std::string, std::string> buildMap(std::string	bcs_in){
  std::map<std::string, std::string>  bc_map;
  std::ifstream	bc_file ;
  bc_file.open(bcs_in) ;
  std::string line ;
  std::vector<std::string> line_elements ;
  // for each line put the first field map value, with the second as key
  while (std::getline(bc_file, line) )
  {
    splitName(line, '\t', line_elements) ;
    bc_map[line_elements[1]] = line_elements[0];
    line_elements.clear() ;
  }
  bc_file.close();
  return bc_map ;
}

int hamming_distance(const std::string& fs, const std::string& ss){
  int hm_distance = 0;

  if((fs.length() == ss.length())){

    for(size_t i = 0; i < fs.length(); i++){
      if(!(fs[i] == ss[i])){
        hm_distance++;
      }
    }

  } else {
    // return -1 if strings not same size
    hm_distance = -1;
  }
  return hm_distance;
}

int main(int argc, char *argv[])
{
  if (argc == 1) {
    std::cerr << "Usage: %s <in.bam> <in.bcs> <out.bam \n" << argv[0] ;
    return 1;
  } else if  (argc == 3) {
    std::cerr << "input and output file names required" << std::endl ;
    return 1;
  }

  //generate map of cell bcs and cell names
  std::string bcs = argv[2] ;
  auto bc_map = buildMap(bcs) ;

  //read in bam
  BamReader bfile(argv[1]) ;  //handle opening bam
  char *fn_out = argv[3] ; //output handle
  const bam_hdr_t *header = sam_hdr_read(bfile.in) ; // get header
  bam1_t *aln = bam_init1() ; // initialize empty alignment container
  samFile *fp_out = sam_open(fn_out, "wb") ; //initialize output bam

  sam_hdr_write(fp_out, header) ;
  while (bam_read1(bfile.bz, aln) > 0) { // negative return values are errors

    std::string id = bam_get_qname(aln) ; //get name

    std::vector<std::string> id_elements ;
    char delim = ':';
    splitName(id, delim, id_elements) ; //split up name field
    auto read_seq = id_elements.back() ; //get seq
    auto cbc = read_seq.substr(0, 11) ;
    auto umi = read_seq.substr(11, 10) ;

    //see if cell barcode is in map, if not take first one 1 hamming dist away
    // otherwise report unmatched
    auto bc_it = bc_map.find(cbc);
    std::string bc ;
    if ( bc_it != bc_map.end()) {
      bc = bc_it->second ;
    } else {
      int min_hamming = 1 ;
      for ( auto& bc_seq : bc_map ){
        auto hd = hamming_distance(cbc, bc_seq.first) ;
        if (hd == min_hamming) {
          bc = bc_seq.second ;
          std::cout << cbc << " " << bc_seq.first << std::endl ;
          }
        }
      if (bc.empty()){
        bc = "Cell_unmatched" ;
      }
    }

    //get original name minus the read seq
    id_elements.erase(id_elements.end() - 1) ;
    auto id_str = join(id_elements, ":") ;

    auto cbc_c = cbc.c_str() ;
    auto bc_c = bc.c_str() ;
    auto umi_c = umi.c_str() ;
    bam_aux_append(aln, "CN", 'Z', strlen(bc_c) + 1, (uint8_t *)(bc_c));
    bam_aux_append(aln, "BX", 'Z', strlen(umi_c) + 1, (uint8_t *)(umi_c));
    bam_aux_append(aln, "CB", 'Z', strlen(cbc_c) + 1, (uint8_t *)(cbc_c));
    sam_write1(fp_out, header, aln) ; //write alignment

  }

  bam_destroy1(aln) ;
  sam_close(fp_out) ;

  return 0;
}
