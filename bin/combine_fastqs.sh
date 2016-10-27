#!/usr/bin/env bash
set -u -e -o pipefail

in_dir=$1

files=$(echo $in_dir*/*.fastq)

for file in $files
  do 
  new_file=$(basename $file) 
  cat $file >> $in_dir$new_file
done 

dirs=$(echo $in_dir*/)
for old_dir in $dirs
 do 
 rm -f $old_dir*
 rmdir $old_dir
done  
