#!/bin/bash
set -euo pipefail

trait=${1:?GWAS trait ID is required}
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

Rscript "${script_dir}/summarize_clumped_variants.R" "${trait}"

tail -n+2 data/GWAS/genes.protein_coding.v39.gtf | \
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$7}' | bedtools sort > protein_coding_genes.bed

awk '{if($1 != "chr23") {print}}' "tmp_${trait}_top.bed" | \
bedtools sort > tmp.bed
bedtools closest -t first -a tmp.bed -b protein_coding_genes.bed | \
awk 'BEGIN{FS="\t"}{print $8,$4}' | sort -gk2 | awk '!a[$1]++' - | \
grep -v -f data/GWAS/closest_genes/HLA.genes - > "data/GWAS/closest_genes/${trait}_genes_closest.txt"
rm -f tmp.bed
