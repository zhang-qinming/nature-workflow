
args   <- commandArgs(trailingOnly = T)  
library(data.table)
df1<-fread( "data/Perturbseq/geneLevel/K562GW/limma_logFC_sum.txt", data.table=F, header=T)
row.names(df1)<-df1[,1]
df1<-df1[,-1]

df2<-fread("data/Perturbseq/geneLevel/K562GW/limma_adjP_sum.txt", data.table=F, header=T)
row.names(df2)<-df2[,1]
df2<-df2[,-1]


corresp<-read.table("data/gencode_v41_gname_gid_ALL_sorted_onlyID", header=F, stringsAsFactor=F)
corresp<-corresp[!duplicated(corresp[,1]),]
corresp<-corresp[!duplicated(corresp[,2]),]
row.names(corresp)<-corresp[,1]

TF_list=dir( "data/TF_ChIP/chip_score")
TF_list=unique(sort( sapply(strsplit(as.character(TF_list), "_"), function(x){x[2]}) ))
summary<-data.frame()

for(TF in TF_list){

if(is.element(TF, colnames(df1))){
df1_f<-df1[,TF,drop=F]
df2_f<-df2[,TF,drop=F]

up_DEG<-intersect(row.names(df1_f)[df1_f[,1]>0], row.names(df2_f)[df2_f[,1]<0.1])
dn_DEG<-intersect(row.names(df1_f)[df1_f[,1]<0], row.names(df2_f)[df2_f[,1]<0.1])

up_DEG=na.omit(corresp[up_DEG, 2])
dn_DEG=na.omit(corresp[dn_DEG, 2])

BG<-na.omit(corresp[row.names(df1_f), 2])


file_list=dir( "data/TF_ChIP/chip_score")
file_list=file_list[grep(paste0("_", TF, "_"), file_list)]
for(FILE in file_list){
ID=sapply(strsplit(as.character(FILE), "_"), function(x){x[1]})
score<-read.table(paste0("data/TF_ChIP/chip_score/", FILE), header=T, stringsAsFactor=F)

for(G in c(1000, 5000, 10000, 50000, 100000, 500000, 1000000 )){
score_f<-score[is.element(score$G, G),]
score_f<-score_f[,c("GENE", "SCORE_qValue")]

##MWU
if(length(up_DEG)>4){
MWU_up<-wilcox.test(score_f[,2][is.element(score_f[,1], up_DEG)], score_f[,2][is.element(score_f[,1], BG[!is.element(BG, up_DEG)])], alternative="g")$p.value
} else {
  MWU_up<- NA  
}

if(length(dn_DEG)>4){
MWU_dn<-wilcox.test(score_f[,2][is.element(score_f[,1], dn_DEG)], score_f[,2][is.element(score_f[,1], BG[!is.element(BG, dn_DEG)])], alternative="g")$p.value
} else {
  MWU_dn<- NA  
}

tmp<-data.frame(chip_ID=ID, TF=TF, SCORE="chip_score", G=G, N_DEG_up=length(up_DEG), N_DEG_dn=length(dn_DEG), MWU_up=MWU_up, MWU_dn=MWU_dn )
summary<-rbind(summary, tmp)
}}


write.table(summary, paste0("data/TF_ChIP/chip_score_KOeffects_MWUtest.txt"), row.names=F, sep="\t", quote=F)
