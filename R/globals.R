# globals shared across markdown docs
library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(openxlsx)
library(Matrix)
library(viridis)

#### Paths ####

project_dir <- path.expand("~/Projects/subsampled_scRNA")
data_dir <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")
docs_dir <- file.path(project_dir, "docs")
db_dir <- file.path(project_dir, "dbases")

# vector of figure paths
figs_dir <-  file.path(results_dir, "Figures") %>%
  dir(pattern = "Figure_[1-4]$",
      include.dirs = TRUE,
      full.names = T)

##### colors #####

# Hexcodes from colorblindr
# palette from  http://jfly.iam.u-tokyo.ac.jp/color/.
palette_okabeito <- c("#E69F00", 
                      "#56B4E9", 
                      "#009E73", 
                      "#F0E442", 
                      "#0072B2", 
                      "#D55E00", 
                      "#CC79A7", 
                      "#999999")
color_palette <- palette_okabeito[5:6]
##### Functions ####

#' When writing out excel workbooks using openxlsx::write.xlsx()
#' this function will set the class attributes for a column, which
#' enforces a column type in the resulting xlsx file. 
#' Useful for avoid gene names being clobbered to dates and 
#' setting scientific number formatting
set_xlsx_class <- function(df, col, xlsx_class){
  for(i in seq_along(col)){
    class(df[[col[i]]]) <- xlsx_class
  }
  df
}

#' write tsv output as gzipped
write_gztsv <- function(df, name, ...){
  if(stringr::str_detect(name, ".gz$")){
    uncompressed_name <- stringr::str_replace(name, ".gz$", "")
  } else {
    uncompressed_name <- name
  }
  readr::write_tsv(df, uncompressed_name, ...)
  system(paste0("gzip -f ", uncompressed_name))
}

#' plot continuous feature on a tSNE
plot_feature <- function(seurat_obj,
                         ident = NULL,
                         gene = NULL,
                         plot_dat = NULL,
                         plot_col = NULL,
                         pt.size = 0.001,
                         pt.alpha = 0.25,
                         label_text = FALSE,
                         label.size = 6,
                         .cols = NULL,
                         legend_name = NULL,
                         cell_filter = NULL,
                         palette_type = "cloupe",
                         col_pal = "Reds",
                         max_y = NULL,
                         ...){
  
  mdata <- seurat_obj@meta.data %>% tibble::rownames_to_column("cell")
  
  tsne_dat <- seurat_obj@dr$tsne@cell.embeddings %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("cell")
  
  tsne_dat <- left_join(mdata, tsne_dat, by = "cell")
  
  if (!is.null(cell_filter)){
    tsne_dat <- dplyr::filter(tsne_dat,
                              cell %in% cell_filter)
  }
  
  if (!is.null(gene)) {
    gene_dat <- FetchData(seurat_obj, gene) %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column("cell")
    tsne_dat <- left_join(tsne_dat, gene_dat, by = "cell")
  }
  
  if (!is.null(plot_dat)){
    tsne_dat <- left_join(tsne_dat, plot_dat, by = "cell")
  }
  
  if (!is.null(ident)){
    color_aes <- ident
  } else if (!is.null(plot_dat)){
    color_aes <- plot_col
  } else {
    color_aes <- gene
  }
  
  p <- ggplot(tsne_dat, 
              aes(tSNE_1, tSNE_2)) +
    geom_point(aes_string(color = color_aes),
               size = pt.size,
               alpha = pt.alpha)
  
  ## handle legend limit 
  
  if (is.null(max_y) & is.null(ident)) {
    max_y <- c(0, max(tsne_dat[[color_aes]]))
  } else if (!is.null(ident)){
    max_y <- c(NA, NA)
  }
  
  ## handle colors
  if (is.null(.cols)){
    if (palette_type == "viridis") {
      p <- p + viridis::scale_color_viridis(discrete = F,
                                   direction = -1,
                                   option = col_pal,
                                   limits = max_y, 
                                   name = legend_name)
    } else if (palette_type == "brewer") {
      p <- p + scale_color_distiller(limits = max_y,
                                     palette = col_pal,
                                     direction = 1, 
                                     name = legend_name)
    } else if (palette_type == "cloupe") {
      cols <- RColorBrewer::brewer.pal(11, "RdGy")[c(1:5, 7)]
      
      p <- p + scale_color_gradientn(limits = max_y,
                                     colors = rev(cols), 
                                     name = legend_name)
    }
  }
  p
}


#' fill a matrix with missing rownames present in 
#' another matrix
#' 
#' @param xmat reference matrix from which to take rownames
#' @param ymat matrix to fill
#' @param fill_value value to fill values in ymat
standardize_rows <- function(xmat,
                             ymat,
                             fill_value = 0L){
  
  stopifnot(all(colnames(xmat) == colnames(ymat)))
  
  ref_rows <- rownames(xmat)
  ref_cols <- colnames(xmat)
  ncols <- length(ref_cols)
  rows_to_add <- setdiff(ref_rows, rownames(ymat))
  zero_mat <- matrix(fill_value, 
                     ncol = ncols, 
                     nrow = length(rows_to_add),
                     dimnames = list(rows_to_add, ref_cols))
  
  res <- rbind(ymat, 
               zero_mat)
  res
}


#' Combine and compare two umi flat files
#'
#'
compare_umis <- function(path_to_ctrl,
                         path_to_test,
                         return_summary = F,
                         strip_10x_suffix = T,
                         cells_exclude = "Cell_unmatched"){
  
  ## umi seqs should be produced by ./get_molecule_info
  ctrl_umi_seqs <- read_tsv(path_to_ctrl,
                            col_names = c("barcode_10x", 
                                          "umi_molecule", 
                                          "count")) %>% 
    filter(!barcode_10x %in% cells_exclude)
  
  test_umi_seqs <- read_tsv(path_to_test,
                            col_names = c("barcode_10x", 
                                          "umi_molecule", 
                                          "count")) %>% 
    filter(!barcode_10x %in% cells_exclude) 
  
  umi_seqs <- full_join(ctrl_umi_seqs, 
                        test_umi_seqs, 
                        by = c("barcode_10x", "umi_molecule"))
  if (strip_10x_suffix) {
    umi_seqs <- mutate(umi_seqs,
                       barcode_10x = str_replace(barcode_10x,
                                                 "-[0-9]$", 
                                                 ""))
  }
  
  if (return_summary) {
    umi_seqs %>% 
      mutate(new_umi = ifelse(is.na(count.x) & !is.na(count.y), 
                              1L, 
                              0L),
             not_detected_umi = ifelse(!is.na(count.x) & is.na(count.y),
                                       1L,
                                       0L),
             shared_umi = ifelse(!is.na(count.x) & !is.na(count.y),
                                 1L,
                                 0L)) %>% 
      group_by(barcode_10x) %>% 
      summarize(new_umis = sum(new_umi),
                not_detected_umis = sum(not_detected_umi),
                shared_umis = sum(shared_umi))
  } else {
    umi_seqs
  }
}

#' simple class to hold info for each experiment
create_sc_obj <- function(umi_df,
                          read_df,
                          cell_mdata_df){
  x <- list()
  class(x) <- "resampled-set"
  x$umis <- umi_df
  x$reads <- read_df
  x$meta_dat <- cell_mdata_df
  return(x)
}


tidy_to_matrix <- function(df){
  df <- as.data.frame(df)
  rownames(df) <- df[, 1]
  df[, 1] <- NULL
  mat <- as.matrix(df)
  mat <- as(mat, "sparseMatrix")   
  return(mat)
}

#' keep both tidy and matrix objs
generate_matrices <- function(sc_obj){
  sc_obj$umi_matrix <- tidy_to_matrix(sc_obj$umis)
  sc_obj$read_matrix <- tidy_to_matrix(sc_obj$reads)
  sc_obj
}

#' normalize by library size (Reads per Million)
norm_libsize <- function(sc_obj){
  sc_obj$norm_umi <- 1e6 * sweep(sc_obj$umi_matrix, 2, 
                                 sum(as.vector(sc_obj$umi_matrix)), "/")
  sc_obj$norm_reads <- 1e6 * sweep(sc_obj$read_matrix, 2, 
                                   sum(as.vector(sc_obj$read_matrix)), "/")
  sc_obj
}

add_metadata <- function(sc_obj, dat, by_col = "barcode_10x"){
  if (is.vector(dat)){
    new_colname <- deparse(substitute(dat))
    df <- data_frame(!!new_colname := dat)
    df[[new_colname]] <- dat
    df[[by_col]] <- names(dat)
    sc_obj$meta_dat <- left_join(sc_obj$meta_dat,
                                 df,
                                 by = by_col)
    
  } else if (is.data.frame(dat)) {
    sc_obj$meta_dat <- left_join(sc_obj$meta_dat,
                                 dat,
                                 by = by_col)
  }
  sc_obj
}

compute_summaries <- function(sc_obj, by_column = "barcode_10x"){
  ## raw counts
  total_umis <- colSums(sc_obj$umi_matrix)
  names(total_umis) <- colnames(sc_obj$umi_matrix)
  total_reads <- colSums(sc_obj$read_matrix)
  names(total_reads) <- colnames(sc_obj$read_matrix)
  
  ## norm counts
  norm_total_umis <- colSums(sc_obj$norm_umi)
  names(norm_total_umis) <- colnames(sc_obj$norm_umi)
  norm_total_reads <- colSums(sc_obj$norm_reads)
  names(norm_total_reads) <- colnames(sc_obj$norm_reads)
  
  sc_obj <- add_metadata(sc_obj, total_umis, by_col = by_column)
  sc_obj <- add_metadata(sc_obj, total_reads, by_col = by_column)
  sc_obj <- add_metadata(sc_obj, norm_total_umis, by_col = by_column)
  sc_obj <- add_metadata(sc_obj, norm_total_reads, by_col = by_column)
  
  ## compute cDNA duplication rate 
  sc_obj$meta_dat$cDNA_duplication <- 1 - (sc_obj$meta_dat$total_umis /
                                             sc_obj$meta_dat$total_reads)
  
  sc_obj
}

umis_to_genes <- function(umipath, 
                          strip_10x_suffix = T, 
                          cells_to_exclude = c("Cell_unmatched")){
  umis <- read_tsv(umipath,
                   col_names = c("barcode_10x", 
                                 "umi_molecule", 
                                 "count")) %>% 
    filter(barcode_10x != cells_to_exclude) %>% 
    mutate(barcode_10x = str_replace(barcode_10x, "-[0-9]$", ""))
  
  mol_fields <- str_count(umis$umi_molecule[1], "::")
  
  if(mol_fields == 2 ){
    umis <- separate(umis, umi_molecule, 
                     into = c("seq", "genome", "gene"),
                     sep = "::") %>% 
      mutate(gene = str_c(genome, "::", gene))
  } else if (mol_fields == 1){
    umis <- separate(umis, umi_molecule, 
                     into = c("seq", "gene"),
                     sep = "::")
  } else {
    stop("separator :: missing from umi_molecule field")
  }
  
  reads <- select(umis, 
                  barcode_10x, 
                  gene,
                  count)
  
  reads <- group_by(reads, 
                    barcode_10x, gene) %>% 
    summarize(counts = sum(count))
  
  reads <- spread(reads, barcode_10x, counts, 
                  fill = 0L)

  reads
}

