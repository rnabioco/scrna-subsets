#' Find identical molecules in a bam file
#' 
#' @description This function uses a cell barcode, umi sequence, overlapping gene, and
#' the genomic position of an alignment to find identical alignments in a bam file. 
#' This is useful to determine what proportion of new UMIs in a subsampling experiment 
#' are identical to molecules detected in other cells in the original library. 
#' Requires a positionally indexed BAM file and bam tags with the UMI seq, and gene id.
#' 
#' @param umi_df data.frame with cell, umi sequence, gene, and position of alignment, 
#' @param bampath path to bam file to find identical molecules

find_molecules <- function(df, bampath,
                           umi_col = "umi",
                           gene_col = "gene",
                           chrom_col = "chrom",
                           pos_col = "pos",
                           umi_tag = "BO",
                           gene_tag = "XT") {
  
  if (!all(c(umi_col, gene_col, pos_col) %in% colnames(df))){
    stop(paste0(umi_col, " or ", gene_col, " or ", pos_col, 
                " not found in input df"))
  }
  
  bampath <- path.expand(bampath)
  if(!file.exists(bampath)){
    stop("bam file not found")
  }
  
  res <- vector(length = nrow(df), mode = "integer")
  
  for(i in seq_along(rownames(df))){
    df_row <- df[i, , drop = T]
    positions <- df_row[[pos_col]][[1]]
    region_strs = str_c(df_row[[chrom_col]],
                       ":",
                       positions)
    
    alignments <- map_dfr(region_strs, 
        ~kentr::bam_to_df(bampath, 
                       tags = c(umi_tag, gene_tag), 
                       region = .x))
    
    ## count number of shared molecules
    molecule_found <- sum(alignments[[gene_tag]] == df_row[[gene_col]] &
                          alignments[[umi_tag]] == df_row[[umi_col]] &
                          alignments[["start"]] %in% positions)
    
    res[i] <- molecule_found
  }
  return(res)
}
