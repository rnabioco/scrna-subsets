#include "subsampling.h"

int usage(){
  std::cerr << "./get_molecule_info <in.bam>" << std::endl ;
  return 1 ;
}

void write_output(std::string current_cell,
                  std::map<std::string, int> umimap ){
  for(auto const &entry : umimap) {
    std::cout << current_cell << "\t"
              << entry.first << "\t"
              << entry.second << std::endl ;
  }
}

void get_molecules(std::string filename){
  // make map to track each unique umi::gene combo
  std::map<std::string, int> umimap ;
  
  BamReader bfile(filename, false) ;
  bam_hdr_t *h = sam_hdr_read(bfile.in) ; // get header
  bam1_t *aln = bam_init1() ;
  
  // get first alignment and set cell group
  int l = bam_read1(bfile.bz, aln) ;
  std::string current_cell;
  if (l > 0){
    Molecule molecule(aln) ;
    current_cell = molecule.cn ;
    umimap[molecule.umikey] += 1 ;
  }  
  
  while (bam_read1(bfile.bz, aln) > 0){
    Molecule molecule(aln) ;
    if(molecule.cn == current_cell){
      umimap[molecule.umikey] += 1 ; 
    } else {
      write_output(current_cell, umimap) ;
      current_cell = molecule.cn ;
      umimap.clear() ;
      umimap[molecule.umikey] += 1 ; 
    }
  }
  // write out final map
  write_output(current_cell, umimap) ;
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
