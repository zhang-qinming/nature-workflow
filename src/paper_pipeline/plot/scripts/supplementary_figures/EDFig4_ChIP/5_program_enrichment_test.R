args   <- commandArgs(trailingOnly = T)   # enable reading factors from linux command line

options("scipen"=10) ##disable 100000 to 1e6
dir.create("data/TF_ChIP/cNMF_enrichment", showWarnings=F)

P = as.character( args[1] ) ##1-60
K = if (length(args) >= 2) as.numeric(args[2]) else 60

TF_list<-read.table("data/TF_ChIP/selected_scores_for_prediction.txt", header=T, stringsAsFactor=F)

##data loading
GEP<-read.table(paste0("data/Perturbseq/cNMF/K562GW/cNMF_all.gene_spectra_score.k_", K, ".dt_0_5.txt"), header=T, stringsAsFactor=F)
GEP<-t(GEP)
colnames(GEP)<-paste0("P", c(1:K))
corresp<-read.table("data/gencode_v41_gname_gid_ALL_sorted_onlyID", header=F, stringsAsFactor=F)
corresp<-corresp[!duplicated(corresp[,1]),]
row.names(corresp)<-corresp[,1]
GEP<-data.frame(GENE=corresp[row.names(GEP), 2], GEP)

for(i in 1:K){
tmp<-read.table(paste0("data/Perturbseq/cNMF_regulation/K562GW/K",  K, "_program", i, "_perturb_effects.txt"), header=T)

tmp1<-tmp[,c(1,2)]
tmp2<-tmp[,c(1,3)]
tmp3<-data.frame(GENE=tmp$GENE, Z=sign(tmp[,2])*qnorm(tmp[,3]/2, lower.tail=F))

colnames(tmp1)<-c("GENE", paste0("P", i))
colnames(tmp2)<-c("GENE", paste0("P", i))
colnames(tmp3)<-c("GENE", paste0("P", i))

if(i==1){
GEP_reg_beta<-tmp1
GEP_reg_P<-tmp2
GEP_reg_Z<-tmp3
} else {
GEP_reg_beta<-merge(GEP_reg_beta, tmp1, by="GENE")
GEP_reg_P<-merge(GEP_reg_P, tmp2, by="GENE")
GEP_reg_Z<-merge(GEP_reg_Z, tmp3, by="GENE")
}
}

##

TOP_N=300
summary<-data.frame()

for(TF in unique(TF_list$TF)){

ID=TF_list[is.element(TF_list$TF, TF), "chip_ID"]
G=TF_list[is.element(TF_list$TF, TF), "G"]
TYPE=TF_list[is.element(TF_list$TF, TF), "TYPE"]

FILE=dir("data/TF_ChIP/chip_score")
FILE=FILE[grep(ID, FILE)]

TFscore<-read.table(paste0("data/TF_ChIP/chip_score/", FILE), header=T, stringsAsFactor=F)
TFscore<-TFscore[is.element(TFscore$G, G),]

TFscore<-TFscore[,c("GENE", "SCORE_qValue")]
colnames(TFscore)<-c("gene", "SCORE")

BG<-intersect(TFscore$gene, GEP$GENE)


reg_Z=GEP_reg_Z[is.element(GEP_reg_Z$GENE, TF), paste0("P", P)]
df1<-GEP[,c("GENE", paste0("P", P))]
colnames(df1)[2]<-"GEP_score"
df1<-df1[order(df1$GEP_score, decreasing=T),]
genes<-df1$GENE[1:TOP_N]

MWU_P<-wilcox.test(TFscore$SCORE[is.element(TFscore$gene, genes)], TFscore$SCORE[is.element(TFscore$gene, BG[!is.element(BG, genes)])], alternative="g")$p.value

hoge<-data.frame(Program=P, TOP_N=TOP_N, TF=TF, TF_type=TYPE, TF_KD_Z=reg_Z, TFscore_program_MWU_P=MWU_P)
summary<-rbind(summary, hoge)
}

summary<-summary[order(summary$TFscore_program_MWU_P),]
write.table(summary, paste0("data/TF_ChIP/cNMF_enrichment/P", P, "_ChIPScoreenrichment.txt"), row.names=F, sep="\t", quote=F)
