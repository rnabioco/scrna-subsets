#!/usr/bin/env bash
set -u -e -o pipefail

USAGE="USAGE: $(basename $0) <FASTQ_FILE> <CHUNK_SIZE> <PREFIX> <SUFFIX>"

if [ $# -ne 4 ]
then
  echo $USAGE
  exit
fi


infile=$1
chunks=$2
prefix=$3
suffix=$4


gunzip -c $infile | \
  split -l $chunks - $prefix"___"

counter=0
for file in $prefix"___"*
  do
  gzip $file
  counter=$((counter + 1))
  mv $file".gz" ${prefix}"_"${counter}${suffix} 
done

