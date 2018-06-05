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


#bampath <- "/Users/kriemo/Projects/subsampled_scRNA/data/lna_pulldown/star/alignments/sorted_large.bam"
#full_df <- read_tsv("data/lna_pulldown/umis/umigroups_positions.txt.gz", 
#               col_names = c("cell_id", "umi", "pos"),
#               col_types = "ccc") 
#
#full_df <- separate(full_df, umi, 
#                    into = c("umi", "chrom", "species", "gene"), sep = "::") %>% 
#  mutate(gene = str_c(species, "::", gene)) %>% 
#  select(cell_id, umi, gene, chrom, pos)
#
#full_df <- mutate(full_df[1:100000, ], 
#             pos = str_replace(pos, 
#                               ",$", "") %>% str_split(., ",") %>% 
#               map(., ~as.integer(.x)))
#
#df <- full_df[1:100, ]
#
#find_molecules(df, bampath)
#
#
#full_df <- read_tsv("data/lna_control/umis/umigroups_positions.txt.gz", 
#                    col_names = c("cell_id", "umi", "pos"),
#                    col_types = "ccc") 
#
#full_df <- separate(full_df, umi, 
#                    into = c("umi", "chrom", "species", "gene"), sep = "::") %>% 
#  mutate(gene = str_c(species, "::", gene)) %>% 
#  select(cell_id, umi, gene, chrom, pos)
#
#full_df <- mutate(full_df[1:100000, ], 
#                  pos = str_replace(pos, 
#                                    ",$", "") %>% str_split(., ",") %>% 
#                    map(., ~as.integer(.x)))
#
#df <- full_df %>% 
#  dplyr::group_by(cell_id) %>% 
#  slice(1:100) %>% 
#  head(400) %>% 
#  ungroup()
#  
#system.time(find_molecules(df, bampath))
#
#
#library(doParallel)
#library(itertools)
#
#n <- 6
#registerDoParallel(n)
#ris <- system.time(foreach(i = split(df, df$cell_id),
#              .packages = c("kentr", "tidyverse"),
#              .export = c("find_molecules", "bampath"))  %dopar% {
#               find_molecules(i, bampath)})
