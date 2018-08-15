#include "subsampling.h"
#include <memory>
extern "C" {
#include "kseq.h"
}


KSEQ_INIT(gzFile, gzread)

// rtrim string
inline std::string& rtrim(std::string& s, std::string t = " \t\n\r\f\v")
{
  s.erase(s.find_last_not_of(t) + 1);
  return s;
}

void splitName(const std::string &name, char delim, std::vector<std::string> &elems) {
  std::stringstream ss;
  ss.str(name);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
}

typedef std::vector<std::shared_ptr<std::ofstream> > fileNames ;

// Make map with barcode as key and file open objects as value
std::map<std::string, fileNames> getBarcodeMap(std::string fn_in, std::string prefix = ""){
  std::map<std::string, fileNames>  bc_map;
    std::ifstream	bc_file ;
    bc_file.open(fn_in) ;
    std::string line ;

    while (std::getline(bc_file, line) )
    { 
      std::ofstream	s ;
      auto bc = rtrim(line) ; 
      
      fileNames fqs ;
      auto r1_fn = prefix + bc + "_R1.fastq" ;
      auto r2_fn = prefix + bc + "_R2.fastq" ;
      
      // fstreams are not copyable, need to work with pointers
      auto r1 = std::make_shared<std::ofstream> (r1_fn); 
      auto r2 = std::make_shared<std::ofstream> (r2_fn) ; 
      
      fqs.push_back(r1) ;
      fqs.push_back(r2) ;
      bc_map[bc] = fqs;
      
    }
    bc_file.close();
    return bc_map ;
}

int main(int argc, char *argv[])
  {
    gzFile fp, fp2;
    kseq_t *seq_1, *seq_2;
    std::string bc_fn ;
    std::string outpre ; 
    if (argc == 1) {
      fprintf(stderr, "Usage: %s <in.seq1> <in.seq2> <barcodes.txt> <output prefix> \n", argv[0]);
      return 1;
    }
    fp = gzopen(argv[1], "r"); 
    fp2 = gzopen(argv[2], "r"); 
    bc_fn = argv[3];

    if (argc == 5){
      outpre = std::string(argv[4]) ;
    } else {
      outpre = std::string("") ;
    }
    auto bc_map = getBarcodeMap(bc_fn, outpre) ;
    
    seq_1 = kseq_init(fp); 
    seq_2 = kseq_init(fp2); 
    int l ; 
    
    while ((l = kseq_read(seq_1)) >= 0) {
      int k = kseq_read(seq_2) ;
      if (k < 0) {
        std::cerr << "incorrect fq2 length" << std::endl;
        return -1;
      }
      
      std::string fq1_name = seq_1->name.s ;
      std::string fq2_name = seq_2->name.s ;
      if (fq1_name != fq2_name) {
        std::cerr << "name mismatch between read_1 " << fq1_name
                  << " and read 2 " << fq2_name << std::endl ;
      }
      
      std::vector<std::string> name_fields ; 
      splitName(fq1_name, '_', name_fields) ;
      
      std::string bc = name_fields[1] ;
      if (bc_map.find(bc) == bc_map.end()) { continue ; }
      auto out_fns = bc_map[bc] ;
      
      *out_fns[0] << "@" << fq1_name << std::endl
                << seq_1->seq.s << std::endl
                << "+" << std::endl
                << seq_1->qual.s << std::endl ;
      
      *out_fns[1] << "@" << fq2_name << std::endl
                 << seq_2->seq.s << std::endl
                 << "+" << std::endl
                 << seq_2->qual.s << std::endl ;
    }
    
    // check if any records remain in seq_2
    int k = kseq_read(seq_2) ;
    if ( k >= 0 ) {
      std::cerr << "incorrect fq1 length" << std::endl;
      return -1;
    }
    
    kseq_destroy(seq_1); 
    kseq_destroy(seq_2); 
    gzclose(fp); 
    gzclose(fp2); 
    
    for(auto &bcs:bc_map) {
      for(auto &fn:bcs.second){
        fn->close() ;
      }
    }
    
    return 0;
  }
