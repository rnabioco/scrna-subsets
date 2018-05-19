#include "subsampling.h"

int usage(){
  std::cerr << "./get_molecule_position <in.bam>" << std::endl ;
  return 1 ;
}

std::string set_to_str(std::set<int> int_set){
  std::ostringstream oss;

  if (!int_set.empty())
  {
    std::copy(int_set.begin(), int_set.end(), std::ostream_iterator<int>(oss, ","));
  } else {
      std::cout << "empty_set" << std::endl ; 
  }

  return oss.str() ;
}


void write_output(std::string cell,
                  std::map<std::string, std::set<int>> posmap){ 
  for(auto const &entry : posmap) {

    std::cout << cell << "\t"
              << entry.first << "\t"
              << set_to_str(entry.second) << std::endl ;
  }
}

//std::string get_chrom(bam1_t *aln,
//                      bam_hdr_t *header){
//  
//    int32_t tid = aln->core.tid;
//    char *chrom = header->target_name[tid] ;
//    return std::string(chrom) ; 
//}

void get_molecules(std::string filename){
  // make map to track each unique umi::gene and their posiions
  std::map<std::string, std::set<int>> posmap ;

  BamReader bfile(filename, false) ;
  bam_hdr_t *h = sam_hdr_read(bfile.in) ; // get header
  bam1_t *aln = bam_init1() ;
  
  // get first alignment and set cell group
  int l = bam_read1(bfile.bz, aln) ;
  std::string current_cell;
  if (l > 0){
    PosMolecule molecule(aln, h) ;
    current_cell = molecule.cn ;
    posmap[molecule.umikey].insert(molecule.pos); 
  }  
  
  while (bam_read1(bfile.bz, aln) > 0){
    PosMolecule molecule(aln, h) ;
    if(molecule.cn == current_cell){
      posmap[molecule.umikey].insert(molecule.pos); 
    } else {
      write_output(current_cell, posmap) ;
      current_cell = molecule.cn ;
      posmap.clear() ;
      posmap[molecule.umikey].insert(molecule.pos); 
    }
  }
  // write out final map
  write_output(current_cell, posmap) ;
}  
  
int main(int argc, char *argv[]){
  if (argc == 1) {
    return usage() ;
  } 
  
  std::string filename(argv[1]) ;
  
  try {
    get_molecules(filename) ;
  } catch (const std::runtime_error& e) {
    std::cerr << e.what() ;
  }
}
