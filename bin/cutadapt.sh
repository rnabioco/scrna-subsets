#!/usr/bin/env bash


dirname=$1

for file in $dirname/*.fastq
do 
cutadapt -a "VNA{100}" \
  -o ${file/.fastq/_trimmed.fastq.gz} \
  $file > ${file/.fastq/.txt} 
done


