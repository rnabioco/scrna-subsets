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
                         legend_names = NULL,
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
                                   limits = max_y)
    } else if (palette_type == "brewer") {
      p <- p + scale_color_distiller(limits = max_y,
                                     palette = col_pal,
                                     direction = 1)
    } else if (palette_type == "cloupe") {
      cols <- RColorBrewer::brewer.pal(11, "RdGy")[c(1:5, 7)]
      
      p <- p + scale_color_gradientn(limits = max_y,
                                     colors = rev(cols))
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
