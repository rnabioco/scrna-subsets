#include "subsampling.h"

Molecule::Molecule(bam1_t *aln){
  cn = std::string(bam_aux2Z(bam_aux_get(aln, "CN"))) ;
  bo = std::string(bam_aux2Z(bam_aux_get(aln, "BO"))) ;
  ug = bam_aux2i(bam_aux_get(aln, "UG")) ;
  xt = std::string(bam_aux2Z(bam_aux_get(aln, "XT"))) ;
  umikey = bo + "::" + xt ;
}

PosMolecule::PosMolecule(bam1_t *aln, bam_hdr_t *header) {
  pos = (aln)->core.pos ;
  cn = std::string(bam_aux2Z(bam_aux_get(aln, "CN"))) ;
  bo = std::string(bam_aux2Z(bam_aux_get(aln, "BO"))) ;
  xt = std::string(bam_aux2Z(bam_aux_get(aln, "XT"))) ;
  chrom = std::string(header->target_name[aln->core.tid]) ;
  umikey = bo + "::" + chrom + "::" + xt ;
  
}
  
