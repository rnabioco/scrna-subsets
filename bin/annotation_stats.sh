#!/usr/bin/env bash
set -u -e -o pipefail

USAGE="USAGE: $(basename $0) <BAM_FILE> <exon_db.bed> <intron_db.bed> <ncrna_db.bed>"

if [ $# -ne 4 ]
then
  echo $USAGE
  exit
fi

bam=$1

exon=$2
intron=$3
ncrna=$4


echo -ne $bam"\t"
echo -ne "exon_count\t"
exons=$(samtools view -bF 4 $bam | bedtools intersect -bed -u -a - -b $exon -nonamecheck | wc -l)
echo  $exons

echo -ne $bam"\t"
echo -ne "intron_count\t"
introns=$(samtools view -bF 4 $bam | bedtools intersect -bed -v -a - -b $exon -nonamecheck | \
  bedtools intersect -u -a - -b $intron -nonamecheck | wc -l)
echo $introns

echo -ne $bam"\t"
echo -ne "ncRNA_count\t"
ncRNAs=$(samtools view -bF 4 $bam | bedtools intersect -bed -v -a - -b $exon -nonamecheck | \
  bedtools intersect -v -a - -b $intron -nonamecheck | \
  bedtools intersect -u -a - -b $ncrna  -nonamecheck | wc -l)
echo $ncRNAs


echo -ne $bam"\t"
echo -ne "intergenic_count\t"
intergenic=$(samtools view -bF 4 $bam | bedtools intersect -bed -v -a - -b $exon -nonamecheck | \
  bedtools intersect -v -a - -b $intron -nonamecheck | \
  bedtools intersect -v -a - -b $ncrna  -nonamecheck | wc -l)
echo $intergenic

