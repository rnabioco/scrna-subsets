#!/usr/bin/env bash
#BSUB -J scRNA
#BSUB -o logs/snakemake_%J.out
#BSUB -e logs/snakemake_%J.err
#BSUB -R "select[mem>4] rusage[mem=4]"

set -o nounset -o pipefail -o errexit -x

args=' -q rna -o {log}.out -e {log}.err -J {params.job_name} -R "{params.memory} " -R "select[hname!=compute11] " -n {threads} '

snakemake --drmaa "$args" --snakefile scRNA.snake --jobs 8 \
  --latency-wait 300 --rerun-incomplete  --configfile config_all.yaml
