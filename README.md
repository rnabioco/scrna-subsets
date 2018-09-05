# scrna-subsets

This repo contains a scRNA-seq data processing pipeline and scripts for analysis for transcriptome resampling libraries. RMarkdown documents are provided to produce figures from the following publication.


>Recovery and analysis of transcriptome subsets from pooled single-cell RNA-seq libraries
Kent A. Riemondy Jr., Monica Ransom, Christopher Alderman, Austin E. Gillen, Rui Fu, Jessica Finlay-Schultz, Gregory Kirkpatrick, Jorge Di Paola, Peter Kabos, Carol A. Sartorius, Jay R. Hesselberth
bioRxiv 408740; doi: https://doi.org/10.1101/408740.


### pipeline

[snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline for processing scRNA-seq fastqs to count matrices, assembling TCR seqs, and producing UMI summary flat files. 

### R

This directory contains a few utility R scripts. The `globals.R` file sets up constant variables and loads commonly used packages.

### results

Collection of RMarkdown documents for producing figures. A snakemake pipeline is provided in the `notebook` subdirectory which runs the rmarkdowns in the correct order. 

### bin 

The snakemake data processing pipeline uses some C++ programs to perform fastq and bam processing steps. These are not necessary to run the analysis code. To compile these programs (requires C++11 compatible compiler, tested on linux and macOS):

```bash
cd bin
make
```
