#include "subsampling.h"

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
    std::cerr << "Usage: %s <in.bam> <in.bcs> <out.bam> <cbc_start> <cbc_end> <umi_start> <umi_end> \n" << argv[0] ;
    return 1;
  } else if  (argc != 8) {
    std::cerr << "input and output file names required" << std::endl ;
    return 1;
  }
  
  int cbc_start, cbc_end, umi_start, umi_end ;
  cbc_start = std::stoi(argv[4]) ; 
  cbc_end   = std::stoi(argv[5]) ;
  umi_start = std::stoi(argv[6]) ;
  umi_end   = std::stoi(argv[7]) ;
  
  
  int cbc_start_idx = cbc_start - 1 ;
  int umi_start_idx = umi_start - 1 ;
  
  int cbc_len = cbc_end - cbc_start_idx ;
  int umi_len = umi_end - umi_start_idx ; 
  
  //generate map of cell bcs and cell names
  std::string bcs = argv[2] ;
  auto bc_map = buildMap(bcs) ;

  //read in bam
  BamReader bfile(argv[1]) ;  //handle opening bam
  char *fn_out = argv[3] ; //output handle
  const bam_hdr_t *header = sam_hdr_read(bfile.in) ; // get header
  bam1_t *aln = bam_init1() ; // initialize empty alignment container
  samFile *fp_out = sam_open(fn_out, "wb") ; //initialize output bam
  
  // track seen barcodes and write to stdout
  std::map<std::string, int> seen_bcs ;

  sam_hdr_write(fp_out, header) ;
  while (bam_read1(bfile.bz, aln) > 0) { // negative return values are errors

    std::string id = bam_get_qname(aln) ; //get name

    std::vector<std::string> id_elements ;
    char delim = ':';
    splitName(id, delim, id_elements) ; //split up name field
    auto read_seq = id_elements.back() ; //get seq
    auto cbc = read_seq.substr(cbc_start_idx, cbc_len) ;
    auto umi = read_seq.substr(umi_start_idx, umi_len) ;

    // count seen barcode
    seen_bcs[cbc] += 1 ; 

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
  
  for(auto const& it:seen_bcs) {
    std::cout << it.first << "\t" << it.second << std::endl ; 
  }

  bam_destroy1(aln) ;
  sam_close(fp_out) ;

  return 0;
}
