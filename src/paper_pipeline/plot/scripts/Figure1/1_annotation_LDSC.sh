#!/bin/bash
set -euo pipefail

dir=${1:?annotation bed name is required}
lddir=${2:?LD reference directory is required}
baselinedir=${3:?baseline reference directory is required}
ldsc_py=${4:?ldsc.py executable is required}

wd=$(pwd)

mkdir -p "ANNOTATION_EUR/${dir}"
cp "bed/${dir}.bed" "ANNOTATION_EUR/${dir}/input.bed"

for chr in $(seq 1 22); do
  cd "${wd}/ANNOTATION_EUR/${dir}"

  awk -v CHR="chr${chr}" '{OFS="\t"; if($1==CHR){gsub("chr","", $1); print $0}}' ./input.bed | bedtools sort > "tmp.${chr}.bed"

  N=$(wc -l "tmp.${chr}.bed" | awk '{print $1}')
  if test "$N" = 0; then
    echo -e "${chr}\t1\t2\t0" > "tmp.${chr}.bed"
  fi

  echo -e "CHR\tBP\tSNP\tCM\tSCORE" > "chr${chr}.annot"
  awk 'BEGIN{OFS="\t"}{print $1,$4,$4,$2,$3}' "${lddir}/1000G.EUR.QC.${chr}.bim" | \
    bedtools intersect -loj -a - -b "tmp.${chr}.bed" | \
    awk 'BEGIN{OFS="\t"}{if($9=="."){$9=0} print $1,$2,$4,$5,$9}' | \
    uniq >> "chr${chr}.annot"

  gzip -f "chr${chr}.annot"

  "${ldsc_py}" --l2 \
    --bfile "${lddir}/1000G.EUR.QC.${chr}" \
    --ld-wind-cm 1 \
    --annot "chr${chr}.annot.gz" \
    --out "chr${chr}" \
    --print-snps "${baselinedir}/${chr}.snp"

  rm -f "tmp.${chr}.bed"
done
