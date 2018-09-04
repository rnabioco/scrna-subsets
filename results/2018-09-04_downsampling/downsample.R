# to be run on a cluster with 55 cores
library(doParallel)
library(dplyr)
library(purrr)
library(Seurat)
library(here)
library(readr)
library(tibble)

project_dir <- here::here()
data_dir <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")
og_sobj <- readRDS(file.path(results_dir, 
                             "2018-05-15_pbmc", 
                             "original_sobj.rds"))
sobj <- readRDS(file.path(results_dir, 
                          "2018-05-15_pbmc", 
                          "rs_v2_sobj.rds"))

og_sobj <- SetAllIdent(og_sobj, "cell_labels")
sobj <- SetAllIdent(sobj, "cell_labels")

resampled_file <- file.exists("downsampled_mk_cluster_markers_resampling.txt")
not_resampled_file <- file.exists("downsampled_mk_cluster_markers_no_resampling.txt")

get_downsampled_markers <- function(sobj, og_sobj, no_cores = 5) {
  mks <- sobj@meta.data[sobj@meta.data$cell_labels == "Megakaryocytes", ]
  rs_mks <- mks[mks$resampled == "original cell", ]
  not_re_mks <- mks[mks$resampled == "not resampled", ]
  not_mks <- sobj@meta.data[sobj@meta.data$cell_labels != "Megakaryocytes", ]
  
  n_mks_to_test <- seq(0, nrow(not_re_mks), by = 1)
  not_re_mk_sampled <- map(n_mks_to_test, 
                           ~sample_n(not_re_mks, .x))
  
  sampled_mks <- map(not_re_mk_sampled, 
                     ~c(rownames(.x), 
                        rownames(rs_mks)))
  
  all_cells_minus_not_sampled_mks <- map(sampled_mks,
                                         ~c(.x, rownames(not_mks)))
  
  
  cl_tmp <- makeCluster(no_cores)  
  registerDoParallel(cl_tmp)  
  subsampled_markers <- foreach(i=all_cells_minus_not_sampled_mks, 
                                .packages =
                                  c("Seurat")) %dopar% {
                                    tmp_dat <- SubsetData(sobj, cells.use = i)
                                    markers <- FindMarkers(tmp_dat, 
                                                           "Megakaryocytes",
                                                           only.pos = T)
                                    markers
                                  }
  stopCluster(cl_tmp)
  
  og_rs_mks <- og_sobj@meta.data[og_sobj@meta.data$resampled, ]
  
  sampled_mks <- map(not_re_mk_sampled, 
                     ~c(rownames(.x), 
                        rownames(og_rs_mks)))
  
  all_cells_minus_not_sampled_mks <- map(sampled_mks,
                                         ~c(.x, rownames(not_mks)))
  cl_tmp <- makeCluster(no_cores)  
  registerDoParallel(cl_tmp)  
  subsampled_markers_og <- foreach(i=all_cells_minus_not_sampled_mks, 
                                   .packages =
                                     c("Seurat")) %dopar% {
                                       tmp_dat <- SubsetData(og_sobj, cells.use = i)
                                       markers <- FindMarkers(tmp_dat, 
                                                              "Megakaryocytes",
                                                              only.pos = T)
                                       markers
                                     }
  stopCluster(cl_tmp)
  
  subsampled_markers <- map(subsampled_markers, 
                            ~tibble::rownames_to_column(.x, "gene")) 
  names(subsampled_markers) <- n_mks_to_test
  subsampled_markers <- bind_rows(subsampled_markers, .id = "n_mks")

  
  subsampled_markers_og <- map(subsampled_markers_og, 
                               ~tibble::rownames_to_column(.x, "gene")) 
  
  names(subsampled_markers_og) <- n_mks_to_test
  subsampled_markers_og <- bind_rows(subsampled_markers_og, .id = "n_mks")
  
  lst(subsampled_markers,
       subsampled_markers_og)
}


cl <- makeCluster(10)  
registerDoParallel(cl)  
all_markers <- foreach(i=1:10, 
                       .packages =
                         c("Seurat", "doParallel", 
                           "purrr", "dplyr")) %dopar% {
                             markers <- get_downsampled_markers(sobj, 
                                                                og_sobj, 
                                                                no_cores = 4)
                           }
stopCluster(cl)

subsampled_markers <- map_dfr(all_markers, ~.x[[1]], .id = "replicate") 
subsampled_markers_og <- map_dfr(all_markers, ~.x[[2]], .id = "replicate") 

write_tsv(subsampled_markers,
          "downsampled_mk_cluster_markers_resampling.txt")
R.utils::gzip("downsampled_mk_cluster_markers_resampling.txt")

write_tsv(subsampled_markers_og,
          "downsampled_mk_cluster_markers_no_resampling.txt")
R.utils::gzip("downsampled_mk_cluster_markers_no_resampling.txt")
