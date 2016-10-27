library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(cowplot)

args <- commandArgs(trailingOnly = TRUE)
dat <- args[1]
dat <- read_tsv(dat)
#remove last line with total reads
dat <- dat[1:(nrow(dat) - 1),]
dat_name <- str_replace(args[1], "_barcode_counts.txt", "_barcode_distibution.pdf")

dat$Barcode <- factor(dat$Barcode, levels = 
                        dat[order(-dat$Count),]$Barcode)
barplot <- ggplot(dat, aes(Barcode, Count)) +
  geom_point(stat = "identity") +
  ggtitle("Demultiplexing") 
save_plot("dat_name", barplot)
