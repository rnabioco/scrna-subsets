#! /usr/bin/env bash
#BSUB -n 55
#BSUB -R " select[mem>200] rusage[mem=200] span[hosts=1]"
#BSUB -J get_markers
#BSUB -q rna
#BSUB -e err.txt
#BSUB -o out.txt

Rscript downsample.R

