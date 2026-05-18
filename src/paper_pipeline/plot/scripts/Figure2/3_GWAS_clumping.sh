#!/bin/bash
set -euo pipefail

TRAIT=${1:?GWAS trait ID is required}
dd=${2:?1000G LD reference directory is required}

mkdir -p clump
echo -e "CHR\tPOS\tID\tP" > clump/${TRAIT}.assoc
gunzip -c "data/GWAS/${TRAIT}_irnt.tsv.gz" | awk '
BEGIN {
    OFS="\t"
}
NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "variant") variant_col = i
        if ($i == "pval") pval_col = i
        if ($i == "low_confidence_variant") low_conf_col = i
    }
    if (!variant_col || !pval_col) {
        print "Required columns variant/pval not found in GWAS file header" > "/dev/stderr"
        exit 2
    }
    next
}
{
    if (low_conf_col) {
        low_conf = tolower($(low_conf_col))
        if (low_conf == "true" || low_conf == "1") {
            next
        }
    }

    p = $(pval_col) + 0
    if ($(pval_col) == "" || $(pval_col) == "NA" || p <= 0 || p > 1) {
        next
    }

    split($(variant_col), best, ":")
    if (length(best) < 4) {
        next
    }

    print best[1], best[2], best[1] ":" best[2] "_" best[3] "_" best[4], $(pval_col)
}
' >> "clump/${TRAIT}.assoc"

cd clump

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


rm -f ${TRAIT}.top
rm -f ${TRAIT}_100kbp.hh
rm -f ${TRAIT}_100kbp.log
