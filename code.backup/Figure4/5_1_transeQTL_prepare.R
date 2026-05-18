
library(data.table)
##downloaded from eQTLGen web page
data<-fread("data/2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt.gz", header=T, data.table=F)

##for set of genes, calculate mean Z score of each SNP.
##use top 100 genes for each program.

K=60
GEP<-read.table(paste0("data/Perturbseq/cNMF/K562GW/cNMF_all.gene_spectra_score.k_", K, ".dt_0_5.txt"), header=T, stringsAsFactor=F)
GEP<-t(GEP)
colnames(GEP)<-paste0("K562GW_P", c(1:K))
summary<-GEP

dir.create("program_trans", showWarnings=F)
for(i in c(1:60)){
genes100<-row.names(summary)[!is.na(summary[,i])][ order(summary[,i][!is.na(summary[,i])], decreasing=T)][1:100]

data_f<-data[is.element(data$Gene, genes100),]
data_f<-data_f[,c(1:9)]
write.table(data_f, paste0("program_trans/P", i, "_top100genes_trans.txt"), row.names=F, sep="\t", quote=F)
}
