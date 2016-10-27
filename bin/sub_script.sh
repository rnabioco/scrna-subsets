#!/usr/bin/env bash
#BSUB -J kallisto_idx
#BSUB -n 2
#BSUB -R "select[mem>10] rusage[mem=10] span[hosts=1]"
#BSUB -o logs/stdout_%J.out
#BSUB -e logs/stderr_%J.err

# module load kallisto/0.42.2.1
# module load gcc

./bin/move.sh
