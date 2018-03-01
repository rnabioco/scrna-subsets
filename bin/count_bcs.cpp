#include "subsampling.h"

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
