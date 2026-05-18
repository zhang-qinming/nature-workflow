#!/bin/bash
set -euo pipefail

K="${1:-60}"
trait_posterior_file="${2:-Backman_2021_86.per_gene_estimates.tsv}"
if [[ "$trait_posterior_file" = /* ]]; then
    trait_posterior_path="$trait_posterior_file"
else
    trait_posterior_path="data/LoF/GeneBayes_posterior/${trait_posterior_file}"
fi

mkdir -p loo_topgene_prediction/gene_list

LANG=en_EN sort  data/gencode_v41_gname_gid_ALL_sorted_onlyID  > loo_topgene_prediction/gencode_v41_gname_gid_forjoin

head -1 "data/Perturbseq/cNMF/K562GW/cNMF_all.gene_spectra_score.k_${K}.dt_0_5.txt" | perl -pe 's/\t/\n/g' |\
LANG=en_EN sort | LANG=en_EN join - loo_topgene_prediction/gencode_v41_gname_gid_forjoin | awk '{print $2}' | LANG=en_EN sort | uniq  > loo_topgene_prediction/GEP_genes.txt

awk '{print $1}' "data/Perturbseq/cNMF_regulation/K562GW/K${K}_program1_perturb_effects.txt" | LANG=en_EN sort | uniq  > loo_topgene_prediction/regulator_genes.txt
cat loo_topgene_prediction/GEP_genes.txt loo_topgene_prediction/regulator_genes.txt |  LANG=en_EN sort | uniq > loo_topgene_prediction/all_genes_perturbseq.txt

##in case of MCH
awk '{print $1}' "${trait_posterior_path}" |\
  LANG=en_EN sort | uniq  |  LANG=en_EN join - loo_topgene_prediction/gencode_v41_gname_gid_forjoin | awk '{print $2}' |\
   LANG=en_EN sort |  LANG=en_EN join - loo_topgene_prediction/all_genes_perturbseq.txt > loo_topgene_prediction/gene_list/allgenes.txt
###

cd loo_topgene_prediction/gene_list
split -d -l 20 allgenes.txt
