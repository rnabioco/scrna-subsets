# Results pipeline

The `results` directory contains various RMarkdown documents that reproduce the figures. These can be run using a `snakemake` pipeline that will handle dependencies between different documents. 

## Processed Data

These RMarkdown documents depend on processed data files. These files can either be generated using the pipeline in the `pipeline` directory, or can downloaded from a zenodo archive. 

The Rmarkdown documents expect the processed data files to be located in a specific directory structure. A R script is provided that will download the processed data (~6Gb), and place it into the correct directory structure.

```bash
Rscript dl_data.R
```

This will generate a `data` directory in the base project directory with all the necessary files for the Rmarkdown documnets to run. 

## Requirements

~12Gb of free space and ~8Gb of RAM.

The following R packages are required:
```R
# CRAN packages
tidyverse
cowplot
openxlsx
Matrix
viridis
here
Seurat
rjson
ggrepel

# Bioconductor
GenomicFeatures
GenomicRanges
GenomicAlignments
Gviz
rtracklayer
ComplexHeatmap

#github packages
kentr #devtools::install_github("kriemo/kentr")
clustifyR #devtools::install_github("NCBI-hackathons/clustifyR")
```

Snakemake will be needed to run the pipeline, and the `blastn` executable from `ncbi-blast` will also need to be available in your PATH.

## Run pipeline

```bash
# preview jobs
snakemake -npr

# run pipeline using 3 cores (may take 2-3 hours to complete) 
snakemake -j 3

```
