#!/bin/bash

mkdir -p result_summary

dir_list=$( ls Partitioning/EUR )

echo -e "TRAIT\tCELL\tEnrichment\tZ" > result_summary/Effect_Z_summary.txt
for DIR in $dir_list; do
DIR2=$( echo $DIR | perl -pe 's/_irnt.sumstats.gz//g' )
file_list=$( ls Partitioning/EUR/${DIR} | grep result )

for FILE in $file_list; do
CELL=$( echo $FILE | perl -pe 's/baseline1.1_//g' | perl -pe 's/_peaks.narrowPeak.results//g' | perl -pe 's/_ATAC_ENCODE_forUse.hg19.results//g')
tail -n+2 Partitioning/EUR/${DIR}/${FILE} | head -1 |\
awk -v CELL=$CELL -v GWAS=$DIR2 'BEGIN{FS="\t"; OFS="\t"}{print GWAS,CELL,$5,$10}' >> result_summary/Effect_Z_summary.txt
done
done
