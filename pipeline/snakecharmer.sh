#!/usr/bin/env bash
#BSUB -J scRNA_snakemake
#BSUB -o logs/snakemake_%J.out
#BSUB -e logs/snakemake_%J.err
#BSUB -R "select[mem>4] rusage[mem=4]"

set -o nounset -o pipefail -o errexit -x


args=' -q normal -o {log}.out -e {log}.err -J {params.job_name} -R "{params.memory}" -n {threads}'

snakemake --drmaa "$args" --snakefile Snakefile --jobs 10 \
  --latency-wait 300 --rerun-incomplete  --configfile config_derivative_libs.yaml

