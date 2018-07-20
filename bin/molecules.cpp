#include "subsampling.h"

Molecule::Molecule(bam1_t *aln){
  uint8_t *tag ;
  if ((tag = bam_aux_get(aln, "CN"))) {
    cn = std::string(bam_aux2Z(tag)) ;
  } else {
    has_tags = false ;
    return ;
  }
  
  if ((tag = bam_aux_get(aln, "BO"))) {
    bo = std::string(bam_aux2Z(tag)) ;
  } else {
    has_tags = false ;
    return ;
  }
  
  if ((tag = bam_aux_get(aln, "UG"))) {
    ug = bam_aux2i(tag) ;
  } else {
    has_tags = false ;
    return ;
  }
  
  if ((tag = bam_aux_get(aln, "XT"))) {
    xt = std::string(bam_aux2Z(tag)) ;
  } else {
    has_tags = false ;
    return ;
  } 
  
  has_tags = true ; 
  umikey = bo + "::" + xt ;
}

PosMolecule::PosMolecule(bam1_t *aln, bam_hdr_t *header) {
  pos = (aln)->core.pos ;
  
  uint8_t *tag ;
  if ((tag = bam_aux_get(aln, "CN"))) {
    cn = std::string(bam_aux2Z(tag)) ;
  } else {
    has_tags = false ;
    return ;
  }
  
  if ((tag = bam_aux_get(aln, "BO"))) {
    bo = std::string(bam_aux2Z(tag)) ;
  } else {
    has_tags = false ;
    return ;
  }
  
  if ((tag = bam_aux_get(aln, "XT"))) {
    xt = std::string(bam_aux2Z(tag)) ;
  } else {
    has_tags = false ;
    return ;
  } 
  
  has_tags = true ; 
  chrom = std::string(header->target_name[aln->core.tid]) ;
  umikey = bo + "::" + chrom + "::" + xt ;
  
}
  
