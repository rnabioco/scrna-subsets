#!/usr/bin/env bash

set -x -o pipefail -u -e 


kallisto index -i data/dbases/kallisto_idx/homo_GRCh38.idx \
  data/dbases/Homo_sapiens.GRCh38.rel79.cdna.all.fa.gz

kallisto index -i data/dbases/kallisto_idx/mus_GRCm38.idx \
  data/dbases/Mus_musculus.GRCm38.rel79.cdna.all.fa.gz
 

