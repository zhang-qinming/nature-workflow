
args <- commandArgs(trailingOnly = TRUE)
K <- if (length(args) >= 1) as.numeric(args[1]) else 60
mch_posterior_path <- if (length(args) >= 2) as.character(args[2]) else "data/LoF/GeneBayes_posterior/Backman_2021_86.per_gene_estimates.tsv"
rdw_posterior_path <- if (length(args) >= 3) as.character(args[3]) else "data/LoF/GeneBayes_posterior/Backman_2021_88.per_gene_estimates.tsv"
irf_posterior_path <- if (length(args) >= 4) as.character(args[4]) else "data/LoF/GeneBayes_posterior/Backman_2021_106.per_gene_estimates.tsv"

corresp<-read.table("data/gencode_v41_gname_gid_ALL_sorted_onlyID", header=F, stringsAsFactor=F)
corresp<-corresp[!duplicated(corresp[,1]),]
row.names(corresp)<-corresp[,1]

##MCH
MCH<-read.table(mch_posterior_path, sep="\t", quote="", header=T, stringsAsFactor=F)
MCH<-data.frame(gene=corresp[MCH$ensg, 2], MCH)

##RDW
RDW<-read.table(rdw_posterior_path, sep="\t", quote="", header=T, stringsAsFactor=F)
RDW<-data.frame(gene=corresp[RDW$ensg, 2], RDW)

##Immature ret fra
IRF<-read.table(irf_posterior_path, sep="\t", quote="", header=T, stringsAsFactor=F)
IRF<-data.frame(gene=corresp[IRF$ensg, 2], IRF)

gene_list=c("CAD", "CALR", "POLE", "SUPT5H", "MED17", "ATR", "GTF3C1", "RBBP6")

MCH_f<-MCH[is.element(MCH$GENE, gene_list), c("GENE", "post_mean")]
colnames(MCH_f)[2]<-"MCH"

RDW_f<-RDW[is.element(RDW$GENE, gene_list), c("GENE", "post_mean")]
colnames(RDW_f)[2]<-"RDW"

IRF_f<-IRF[is.element(IRF$GENE, gene_list), c("GENE", "post_mean")]
colnames(IRF_f)[2]<-"IRF"

df<-merge(MCH_f,RDW_f, by="GENE")
df<-merge(df, IRF_f, by="GENE")
row.names(df)<-df$GENE

df<-df[gene_list,]

write.table(df, "selected_genes_gamma.txt", row.names=F, sep="\t", quote=F)

##confirm the direction of program regulation
##data loading
GEP<-read.table(paste0("data/Perturbseq/cNMF/K562GW/cNMF_all.gene_spectra_score.k_", K, ".dt_0_5.txt"), header=T, stringsAsFactor=F)
GEP<-t(GEP)
colnames(GEP)<-paste0("P", c(1:K))
corresp<-read.table("data/gencode_v41_gname_gid_ALL_sorted_onlyID", header=F, stringsAsFactor=F)
corresp<-corresp[!duplicated(corresp[,1]),]
row.names(corresp)<-corresp[,1]
GEP<-data.frame(GENE=corresp[row.names(GEP), 2], GEP)

GEP_reg_beta<-data.frame()
GEP_reg_P<-data.frame()

for(i in 1:K){
tmp<-read.table(paste0("data/Perturbseq/cNMF_regulation/K562GW/K",  K, "_program", i, "_perturb_effects.txt"), header=T)

tmp1<-tmp[,c(1,2)]
tmp2<-tmp[,c(1,3)]

colnames(tmp1)<-c("GENE", paste0("P", i))
colnames(tmp2)<-c("GENE", paste0("P", i))

if(i==1){
GEP_reg_beta<-tmp1
GEP_reg_P<-tmp2
} else {
GEP_reg_beta<-merge(GEP_reg_beta, tmp1, by="GENE")
GEP_reg_P<-merge(GEP_reg_P, tmp2, by="GENE")
}
}


##
GEP_reg_FDR<-data.frame(GENE=GEP_reg_P$GENE, P40=p.adjust(GEP_reg_P[,"P40"], method="BH"), P16=p.adjust(GEP_reg_P[,"P16"], method="BH"),P43=p.adjust(GEP_reg_P[,"P43"], method="BH"), P25=p.adjust(GEP_reg_P[,"P25"], method="BH"))
GEP_reg_FDR_f<-GEP_reg_FDR[is.element(GEP_reg_FDR$GENE, gene_list),]
GEP_reg_beta_f<-GEP_reg_beta[is.element(GEP_reg_beta$GENE, gene_list),c("GENE", "P40", "P16", "P43", "P25")]

write.table(GEP_reg_FDR_f, "selected_genes_regulation_FDR.txt", row.names=F, sep="\t", quote=F)
write.table(GEP_reg_beta_f, "selected_genes_regulation_beta.txt", row.names=F, sep="\t", quote=F)

##draw Fig6B
