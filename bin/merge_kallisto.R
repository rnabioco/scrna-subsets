library(readr)
library(dplyr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
# set  input variable with the kallisto output dirs"

kallisto_dir <- args[1]
message(paste("working on ", kallisto_dir, sep = ""))
files <- list.files(path = kallisto_dir, 
  recursive = T, 
  pattern = "abundance.tsv",
  include.dirs = T,
  full.names = T)

##compact string representation of column types ("ciddd")
df_counts_list <- lapply(files,  readr::read_tsv, col_types = "c--d-")
df_tpm_list <- lapply(files,  readr::read_tsv, col_types = "c---d")

df_names <- stringr::str_split(files, "/")
dir_length <- length(df_names[[1]])
df_names <- lapply(seq(length(df_names)), function(x) df_names[[x]][dir_length - 1])
df_names <- unlist(df_names)

df_names <- c("target_id", df_names)

count_dat <- Reduce(function(...) dplyr::full_join(..., by = "target_id"), df_counts_list)
tpm_dat <- Reduce(function(...) dplyr::full_join(..., by = "target_id"), df_tpm_list)

colnames(count_dat) <- df_names
colnames(tpm_dat) <- df_names

saveRDS(count_dat, paste(kallisto_dir, "kallisto_count_matrix.rds", sep =
""))
saveRDS(tpm_dat, paste(kallisto_dir, "kallisto_tpm_matrix.rds", sep = ""))
