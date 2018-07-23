#include "subsampling.h"

extern "C" {
#include "kseq.h"
}


KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
  gzFile fp, fp2;
  kseq_t *seq_1;
  kseq_t *seq_2;
  int l ;
  int cbc_start, cbc_end, umi_start, umi_end ;
  std::string delim ; 
  
  if (argc < 6) {
    fprintf(stderr, "Usage: %s <in.seq1> <in.seq2> <cbc_start> <cbc_end> <umi_start> <umi_end> <delim> \n", argv[0]);
    return 1;
  }
  fp = gzopen(argv[1], "r"); 
  fp2 = gzopen(argv[2], "r");
  
  cbc_start = std::stoi(argv[3]) ; 
  cbc_end   = std::stoi(argv[4]) ;
  umi_start = std::stoi(argv[5]) ;
  umi_end   = std::stoi(argv[6]) ;
  
  int cbc_start_idx = cbc_start - 1 ;
  int umi_start_idx = umi_start - 1 ;
  
  int cbc_len = cbc_end - cbc_start_idx ;
  int umi_len = umi_end - umi_start_idx ; 
  
  if (argc == 8){
    delim = argv[7] ;
  } else {
    delim = "_" ;
  }
  
  // track seen barcodes and write to stdout
  std::map<std::string, int> seen_bcs ;

  seq_1 = kseq_init(fp); 
  seq_2 = kseq_init(fp2);

  while ((l = kseq_read(seq_1)) >= 0) {

    std::string fq1_seq = seq_1->seq.s ;
    int k = kseq_read(seq_2) ;
    if (k < 0) {
      std::cerr << "incorrect fq2 length" << std::endl;
      return 1;
    }

    std::string fq1_name = seq_1->name.s ;
    std::string fq2_name = seq_2->name.s ;
    if (fq1_name != fq2_name) {
      std::cerr << "name mismatch between read_1 " << fq1_name
                << " and read 2 " << fq2_name << std::endl ;
    }
    
    auto cbc = fq1_seq.substr(cbc_start_idx, cbc_len) ;
    auto umi = fq1_seq.substr(umi_start_idx, umi_len) ;
    std::string  r1_seq_out = cbc + delim + umi ;
    
    std::cout << "@" << fq2_name << delim << r1_seq_out << std::endl
              << seq_2->seq.s << std::endl
              << "+" << std::endl
              << seq_2->qual.s << std::endl ;
    
    seen_bcs[cbc] += 1 ; 
  }

  // check if any records remain in seq_2
  int k = kseq_read(seq_2) ;
  if ( k >= 0 ) {
    std::cerr << "incorrect fq1 length" << std::endl;
    return 1;
  }

  for(auto const& it:seen_bcs) {
    std::cerr << it.first << "\t" << it.second << std::endl ; 
  }
  
  kseq_destroy(seq_1); 
  kseq_destroy(seq_2); 
  gzclose(fp); 
  gzclose(fp2); 

  return 0;
}
