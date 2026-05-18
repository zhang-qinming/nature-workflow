#!/bin/bash

TRAIT=$1 ##loop over traits.

mkdir -p clump
echo -e "CHR\tPOS\tID\tP" > clump/${TRAIT}.assoc
gunzip -c  data/GWAS/${TRAIT}_irnt.tsv.gz | tail -n+2 | awk 'BEGIN{OFS="\t"}{split($1,BEST, ":"); print BEST[1],BEST[2],BEST[1]":"BEST[2]"_"BEST[3]"_"BEST[4], $11}' >> clump/${TRAIT}.assoc

cd clump
module load biology plink/1.90b5.3

dd="path to 1000G LD reference"

plink --bfile ${dd}/all_phase3_onlyEUR  --allow-extra-chr \
--clump  ${TRAIT}.assoc  \
--clump-p1 5e-8 --clump-p2 5e-8 --clump-r2 0.01 --clump-kb 10000 --clump-field P --clump-snp-field ID --clump-verbose \
--out ${TRAIT}

##merge snps within 100kbp

echo -e "CHR\tPOS\tID\tP" > ${TRAIT}.top
grep -A 1 CHR  ${TRAIT}.clumped | grep -v CHR |\
awk 'BEGIN{OFS="\t"}{if(NF==11){print $1,$4,$3,$5}}' | sort -k1,1 -k2,2n >> ${TRAIT}.top

plink --bfile ${dd}/all_phase3_onlyEUR  --allow-extra-chr \
--clump  ${TRAIT}.top  \
--clump-p1 5e-8 --clump-p2 5e-8 --clump-r2 0 --clump-kb 100 --clump-field P --clump-snp-field ID --clump-verbose \
--out ${TRAIT}_100kbp

echo -e "CHR\tPOS\tID\tP" > ${TRAIT}.top.100kbpmerged
grep -A 1 CHR  ${TRAIT}_100kbp.clumped | grep -v CHR |\
awk 'BEGIN{OFS="\t"}{if(NF==11){print $1,$4,$3,$5}}' | sort -k1,1 -k2,2n >> ${TRAIT}.top.100kbpmerged


rm ${TRAIT}.top
rm ${TRAIT}_100kbp.hh
rm ${TRAIT}_100kbp.log
