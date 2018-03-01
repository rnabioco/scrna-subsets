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
  int l;
  if (argc == 1) {
    fprintf(stderr, "Usage: %s <in.seq1> <in.seq2> \n", argv[0]);
    return 1;
  }
  fp = gzopen(argv[1], "r"); 
  fp2 = gzopen(argv[2], "r");
  seq_1 = kseq_init(fp); 
  seq_2 = kseq_init(fp2);

  while ((l = kseq_read(seq_1)) >= 0) {

    std::string fq1_seq = seq_1->seq.s ;
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

    std::cout << "@" << fq2_name << ":" << fq1_seq << std::endl
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

  return 0;
}
