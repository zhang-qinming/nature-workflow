#!/bin/bash

module load biology bedtools
lddir="path to LD reference"
baselinedir="path to baseline reference"
LDSC_path="path to LDSC software"


wd=$( pwd )
dir=$1 ##loop over bed files

mkdir -p ANNOTATION_EUR/${dir}
cp  bed/${dir}.bed ANNOTATION_EUR/${dir}/input.bed


for chr in `seq 1 22`; do

cd ANNOTATION_EUR/${dir}
   
   awk -v CHR=chr$chr '{OFS="\t";
      if($1==CHR){gsub("chr","", $1); print $0}}' ./input.bed | bedtools sort > tmp.$chr.bed
   
   N=$( wc -l tmp.$chr.bed | awk '{print $1}' )
   if test $N = 0 ;then
      echo -e "$chr\t1\t2\t0" > tmp.$chr.bed 
   fi
   
echo -e "CHR\tBP\tSNP\tCM\tSCORE" > chr${chr}.annot
awk 'BEGIN{OFS="\t"}{print $1,$4,$4,$2,$3}'  ${lddir}/1000G.EUR.QC.$chr.bim |\
bedtools intersect -loj -a - -b tmp.$chr.bed | awk 'BEGIN{OFS="\t"}{if($9=="."){$9=0} print $1,$2,$4,$5,$9}' |\
uniq  >> chr${chr}.annot

gzip -f chr${chr}.annot

 ${LDSC_path}/ldsc.py  --l2 \
      --bfile ${lddir}/1000G.EUR.QC.$chr  \
      --ld-wind-cm   1         \
      --annot chr$chr.annot.gz \
      --out   chr$chr          \
      --print-snps ${baselinedir}/$chr.snp
   
   rm -f tmp.$chr.bed
done
