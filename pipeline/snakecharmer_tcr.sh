#!/usr/bin/env bash
#BSUB -J tcrscRNA
#BSUB -o logs/snakemake_%J.out
#BSUB -e logs/snakemake_%J.err
#BSUB -R "select[mem>4] rusage[mem=4]"

set -o nounset -o pipefail -o errexit -x

args=' -q rna -o {log}.out -e {log}.err -J {params.job_name} -R " {params.memory} span[hosts=1] " -n {threads} '


#### load necessary programs ####

# If programs are not all in the path then modify code to load 
# the necessary programs

# load modules
. /usr/share/Modules/init/bash
module load modules modules-init modules-python

module load samtools/1.5
module load bowtie2/2.3.2

#### execute snakemake ####

snakemake --drmaa "$args" \
    --snakefile tcr.snake \
    --jobs 2 \
    --resources all_threads=20 \
    --latency-wait 50 \
    --rerun-incomplete  
