args   <- commandArgs(trailingOnly = T) 

options("scipen"=10) ##disable 100000 to 1e6

k<-as.numeric( args[1] ) ##60
i <-as.numeric( args[2] ) ##1-60
i=i + 1 ##adjust for the first column

metadata=read.csv("data/Perturbseq/metadata/gwps_metadata.csv", row.names = 1)

library(data.table)
count<-fread(paste0("data/Perturbseq/cNMF/K562GW/cNMF_all.usages.k_", k, ".dt_0_5.consensus.txt"), header=T, data.table=F)

GENES<-unique(sort(metadata$gene))
GENES<-GENES[!is.element(GENES, "non-targeting")]

library(dplyr)

control<- row.names( metadata %>% filter(gene %in% "non-targeting") )

N=ncol(count)

summary<-data.frame()

B<- count[is.element(count$cell_barcode, control),c(1,i)] 
dfB<-cbind(B, metadata[B$cell_barcode, c("gem_group", "n_genes", "mitopercent")])
dfB<-data.frame(dfB, GROUP="control")

A<- count[,c(1,i)] 
tmp<-intersect(unique(A$cell_barcode), row.names(metadata))
row.names(A)<-A$cell_barcode
A<-A[tmp,]
A<-cbind(A, metadata[ tmp, c("gem_group", "n_genes", "mitopercent")])



for(GENE in GENES) {
target<- row.names( metadata %>% filter(gene %in% GENE ) )
if(length(target)>10){

dfA<-A[is.element(A$cell_barcode, target),]
dfA<-data.frame(dfA, GROUP="perturb")

if(length(which(is.element(dfB$gem_group, unique(dfA$gem_group))))>0){

df<-rbind(dfA,dfB)
df$gem_group <-factor(as.character(df$gem_group))

colnames(df)[2]<-"exp"
df$exp<-scale(df$exp)

m.lm <- lm(exp ~ GROUP  + gem_group + n_genes + mitopercent, data=df)
lm_es<-summary(m.lm)$coefficients[2,1]
lm_p<-summary(m.lm)$coefficients[2,4]
hoge<-data.frame(GENE=GENE,  lm_es=lm_es, lm_p=lm_p)
summary<-rbind(summary, hoge)

}}}

dir.create(paste0("data/Perturbseq/cNMF_regulation/K562GW"), showWarnings=F)
write.table(summary, paste0("data/Perturbseq/cNMF_regulation/K562GW/K60_program", colnames(count)[i],  "_perturb_effects.txt"), row.names=F, sep="\t", quote=F, append=F)
