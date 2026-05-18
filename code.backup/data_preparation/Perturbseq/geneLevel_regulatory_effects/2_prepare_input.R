
args   <- commandArgs(trailingOnly = T)   

options("scipen"=10) ##disable 100000 to 1e6

CHUNK=as.numeric( args[1] ) ##1,2,3

dir.create("limmatrend", showWarnings=F)
dir.create("limmatrend/tmp_limma_input/", showWarnings=F)

if(CHUNK==1){
range=c(1:60)
}
if(CHUNK==2){
range=c(61:120)
}
if(CHUNK==3){
range=c(121:198)
}

metadata=read.csv("data/Perturbseq/metadata/gwps_metadata.csv", row.names = 1)

library(data.table)
library(dplyr)

count<-fread("data/Perturbseq/K562gwps_raw_count.csv", sep=",", header=T )
control<- row.names( metadata %>% filter(gene %in% "non-targeting") )

if(CHUNK==1){
B<- count[is.element(count$cell_barcode, control),] 
dfB<-cbind(B, metadata[B$cell_barcode, c("gem_group", "n_genes", "mitopercent", "leiden", "gene")])
dfB<-data.frame(dfB, GROUP="control")
write.table(dfB, paste0("limmatrend/tmp_limma_input/control.input"), row.names=F, sep="\t", quote=F)
}

for(AA in range){

library(dplyr)
library(kSamples)


START=50*(AA - 1) + 1 
END=50*AA

if(AA==198){
END=9866
}


GENES<-unique(sort(metadata$gene))
GENES<-GENES[!is.element(GENES, "non-targeting")]
GENES<-GENES[c(START:END)]


target<- row.names( metadata [is.element(metadata$gene, GENES),]) 

A<- count[is.element(count$cell_barcode, target),] 

N=ncol(A)

dfA<-cbind(A, metadata[A$cell_barcode, c("gem_group", "n_genes", "mitopercent", "leiden", "gene")])
dfA<-data.frame(dfA, GROUP="perturb")

write.table(dfA, paste0("limmatrend/tmp_limma_input/", AA, ".input"), row.names=F, sep="\t", quote=F)
}
