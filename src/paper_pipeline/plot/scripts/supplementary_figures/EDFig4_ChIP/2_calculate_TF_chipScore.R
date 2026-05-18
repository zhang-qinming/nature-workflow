args   <- commandArgs(trailingOnly = T)  

chip <- as.character( args[1] ) ##file name for chip-seq narrow peak file
chip2<- gsub(".bed", "", chip)

gene_pos<-read.table("data/genes.protein_coding.v39.gtf", header=T, stringsAsFactor=F)
gene_pos<-gene_pos[,c("chr", "gene", "tss")]
chip_data<-read.table(paste0("data/TF_ChIP/hg19/",  chip), header=F, stringsAsFactor=F)

result<-data.frame()
library(plyr)
library(dplyr)

for(chr in c(1:22)){
tmp_gene_pos<-gene_pos[gene_pos[,1]==paste0("chr",chr),]
df<-merge(tmp_gene_pos, chip_data, by=1)

for(G in c(1000, 5000, 10000, 50000, 100000, 500000, 1000000 )){

tmp<-df[,c(2, 4,5,6,7, 3)]
colnames(tmp)<-c("GENE", "POS1", "POS2", "VALUE1", "VALUE2", "POS3")
tmp$POS3<-as.numeric(as.character(tmp$POS3))

summary <- ddply(.data=tmp, 
                 .(GENE), 
                 summarize, 
                 SCORE_signalValue=sum( VALUE1 * ( exp(-abs( ((POS1+POS2)/2) -POS3 ) /G) )),
                 SCORE_qValue=sum( VALUE2 * ( exp(-abs( ((POS1+POS2)/2) -POS3 ) /G) ))) 
##
summary<-summary[!duplicated(summary$GENE),]

hoge<-data.frame(CHR=paste0("chr", chr), G=G, summary)
result<-rbind(result, hoge)
}
}

write.table(result, paste0("data/TF_ChIP/chip_score/", chip2, "_scores.txt"), row.names=F, sep="\t", quote=F)
