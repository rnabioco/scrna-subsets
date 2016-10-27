library(dplyr)
library(readr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
dat <- read_tsv(input)

#check that row/col descriptors are correct

sanity_check <- nrow(dat) == nrow(unique(dat[ , c(1, 2)]))
if (!sanity_check){
  stop("duplicate row/col id's detected")
}

dat$mod.sample <- str_replace_all(dat$Sample, "_","")
dat$mod.sample <- str_replace_all(dat$mod.sample," ", "")
dat$mod.sample <- str_replace_all(dat$mod.sample,"-", "") 

# extract fastq ids for each cell
dat$fastq_id <- str_c("Cell_", 
      dat$mod.sample, 
      "_", dat$Row, "_", dat$Col, sep = "")
  
dat_min <- dat %>%
  select(fastq_id, Barcode)

sanity_check <- nrow(dat_min) == nrow(dat)
if (!sanity_check){
  stop("error generating fastq_ids")
}

write_tsv(dat_min, args[2], col_names = F)
