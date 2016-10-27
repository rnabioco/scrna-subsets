#!/usr/bin/env bash
set -u -e -o pipefail

USAGE="USAGE: $(basename $0) <List of Fastqs>"

if (( $# == 0 ))
then
  echo $USAGE
  exit
fi



for file in "$@";
do 
  if [[ $file == *.fastq.gz ]]
  then
    count=$(zgrep -c "^@HISEQ:" $file)
    echo -e $file'\t'$count
  else 
    if [[ $file == *.bam ]]
    then 
    count=$(samtools view $file | cut -f 1 | sort | uniq | wc -l)
    echo -e $file'\t'$count
    fi
  fi
done

