#!/usr/bin/env bash
#BSUB -J scRNA
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

module load ucsc/v308 
module load fastqc/0.11.5
module load samtools/1.5
module load STAR/2.5.2a
module load subread/1.6.2

#### execute snakemake ####

snakemake --drmaa "$args" \
    --snakefile 5p.snake \
    --jobs 36 \
    --resources all_threads=36 \
    --latency-wait 50 \
    --rerun-incomplete  \
    --configfile config_human_5p.yaml 

snakemake --drmaa "$args" \
    --snakefile Snakefile \
    --jobs 36 \
    --resources all_threads=36 \
    --latency-wait 50 \
    --rerun-incomplete  \
    --configfile config_human.yaml 

snakemake --drmaa "$args" \
    --snakefile Snakefile \
    --jobs 36 \
    --resources all_threads=36 \
    --latency-wait 50 \
    --rerun-incomplete  \
    --configfile config_mouse_human.yaml 

snakemake --drmaa "$args" \
    --snakefile tcr.snake \
    --jobs 36 \
    --resources all_threads=36 \
    --latency-wait 50 \
    --rerun-incomplete  \
    --configfile config_tcr.yaml 
