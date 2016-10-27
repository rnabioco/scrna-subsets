#!/usr/bin/env bash

set -u -o pipefail -e

input_dir=$1
output=$2

for file in $input_dir/Cell*.fastq 
  do echo $file
  awk -v x=$file 'NR % 4 ==1 {print $1, x}' $file \
  >> $output
done
