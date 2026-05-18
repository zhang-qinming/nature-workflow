#!/bin/bash

mkdir -p data/TF_ChIP
cd data/TF_ChIP

##download K562 ChIP-seq list file from ENCODE ##files.txt

head -1 files.txt > files.txt2
grep bed.gz files.txt >> files.txt2

xargs -L 1 curl -O -J -L < files.txt2

##simplify metadata

file_list=$( ls | grep bed.gz )
rm -f metadata.sum
for FILE in $file_list; do
ID=$( echo $FILE | perl -pe 's/.bed.gz//g' )
grep -w $ID metadata.tsv | grep "bed narrowPeak" | awk -v FILE=$FILE 'BEGIN{FS="\t"; OFS="\t"}{split($23, BEST, "-"); print FILE, $6,BEST[1]}' >> metadata.sum
done

##liftover and rename the file
##column: CHR, START, END, signalValue,qvalue
module load biology bedtools/2.30.0 

mkdir -p hg19
file_list=$( ls | grep bed.gz )
for FILE in $file_list; do
cord=$( grep $FILE  metadata.sum | awk '{print $2}')
TF=$( grep $FILE  metadata.sum | awk '{print $3}')
ID=$( echo $FILE | perl -pe 's/.bed.gz//g' )
FILE2="${ID}_${TF}.bed"

if [ $cord == "hg19" ]; then
gunzip -c $FILE | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$7,$9}' | bedtools sort > hg19/${FILE2}
fi
if [ $cord == "GRCh38" ]; then 
gunzip -c ${FILE} | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$7,$9}' >tmp.bed
 ~/software/liftOver tmp.bed   ~/software/hg38ToHg19.over.chain.gz  hg19/${FILE2} unlifted.bed
rm tmp.bed
bedtools sort -i hg19/${FILE2}  > hg19/${FILE2}_2
mv hg19/${FILE2}_2  hg19/${FILE2}
fi
done

mkdir -p  chip_score
