#!/bin/bash

mkdir -p loo_topgene_prediction/gene_list

LANG=en_EN sort  data/gencode_v41_gname_gid_ALL_sorted_onlyID  > loo_topgene_prediction/gencode_v41_gname_gid_forjoin

head -1 data/Perturbseq/cNMF/K562GW/cNMF_all.gene_spectra_score.k_60.dt_0_5.txt | perl -pe 's/\t/\n/g' |\
LANG=en_EN sort | LANG=en_EN join - loo_topgene_prediction/gencode_v41_gname_gid_forjoin | awk '{print $2}' | LANG=en_EN sort | uniq  > loo_topgene_prediction/GEP_genes.txt

awk '{print $1}' data/Perturbseq/cNMF_regulation/K562GW/K60_program1_perturb_effects.txt | LANG=en_EN sort | uniq  > loo_topgene_prediction/regulator_genes.txt
cat loo_topgene_prediction/GEP_genes.txt loo_topgene_prediction/regulator_genes.txt |  LANG=en_EN sort | uniq > loo_topgene_prediction/all_genes_perturbseq.txt

##in case of MCH
awk '{print $1}'  data/LoF/GeneBayes_posterior/Backman_2021_86.per_gene_estimates.tsv |\
  LANG=en_EN sort | uniq  |  LANG=en_EN join - loo_topgene_prediction/gencode_v41_gname_gid_forjoin | awk '{print $2}' |\
   LANG=en_EN sort |  LANG=en_EN join - loo_topgene_prediction/all_genes_perturbseq.txt > loo_topgene_prediction/gene_list/allgenes.txt
###

cd loo_topgene_prediction/gene_list
split -d -l 20 allgenes.txt
