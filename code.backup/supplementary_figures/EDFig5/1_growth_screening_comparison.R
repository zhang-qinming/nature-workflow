args   <- commandArgs(trailingOnly = T)   # enable reading factors from linux command line
dir.create("data/fitness_compare/", showWarnings=F)

P_target=as.character(args[1]) ##"P1"-"P60"

library(data.table)
df<-read.csv("data/Growth_screening_2014.csv", header=T)  ##adopted from Gilbert et al. 2014.
df<-df[,c(1,2)]
colnames(df)[1]<-"gene"
df<-df[-1,]

shet<-read.table("data/shet_10bins.txt", header=T, stringsAsFactor=F)
shet<-shet[,c("GENE", "shet", "shet_BIN")]
colnames(shet)[1]<-"gene"

K=60

##data loading
GEP<-read.table(paste0("data/Perturbseq/cNMF/K562GW/cNMF_all.gene_spectra_score.k_", K, ".dt_0_5.txt"), header=T, stringsAsFactor=F)
GEP<-t(GEP)
colnames(GEP)<-paste0("P", c(1:K))
corresp<-read.table("data/gencode_v41_gname_gid_ALL_sorted_onlyID", header=F, stringsAsFactor=F)
corresp<-corresp[!duplicated(corresp[,1]),]
row.names(corresp)<-corresp[,1]
GEP<-data.frame(GENE=corresp[row.names(GEP), 2], GEP)

library(dplyr)

program<-GEP[order(GEP[,P_target], decreasing=T), "GENE"][1:100]

df_pro<-df %>% filter(is.element(gene, program))
df_ctrl<-df %>% filter(!is.element(gene, c(program)))

##compare program genes with control genes while matching shet

shet_tmp<-shet[is.element(shet$gene, df_pro$gene),]
shet_tmp_ctrl<-shet[is.element(shet$gene, df_ctrl$gene),]

A<-table(shet_tmp$shet_BIN)

summary<-data.frame(SEED="true", meangamma=mean(as.numeric(df3_pro$CRISPRi)))

for(I in 1:100000){
set.seed(I)
genes<-c()
for(p in 1:length(A)){
genes<-c(genes, sample(shet_tmp_ctrl$gene[is.element(shet_tmp_ctrl$shet_BIN, names(A)[p])], A[p]))
}
hoge<-data.frame(SEED=I, meangamma=mean(as.numeric(df3_ctrl[is.element(df3_ctrl$gene, genes),"CRISPRi"]), na.omit=T))
summary<-rbind(summary, hoge)
}

write.table(summary, paste0("data/fitness_compare/", P_target, "_program_Gamma_distribution.txt"), row.names=F, sep="\t", quote=F)

##from here plot the enrichment vs p value