#!/usr/bin/env bash
#BSUB -J scRNA
#BSUB -o logs/snakemake_%J.out
#BSUB -e logs/snakemake_%J.err
#BSUB -R "select[mem>4] rusage[mem=4]"

set -o nounset -o pipefail -o errexit -x


args=' -q rna -o {log}.out -e {log}.err -J {params.job_name} -R "{params.memory} " -n {threads} '

snakemake --drmaa "$args" --snakefile Snakefile --jobs 20 \
  --latency-wait 300 --rerun-incomplete  --configfile config_phospho_libs.yaml
