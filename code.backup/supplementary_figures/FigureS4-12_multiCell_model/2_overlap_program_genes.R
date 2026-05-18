cell_list=dir("data/Perturbseq/cNMF_regulation")
K=60
dir.create("multiCell_model", showWarnings=F)

summary<-data.frame()
for(CELL in cell_list){
FILE=dir(paste0("data/Perturbseq/cNMF/", CELL))
FILE=FILE[grep("gene_spectra_score", FILE)]

GEP<-read.table(paste0("data/Perturbseq/cNMF/", CELL, "/", FILE), header=T, stringsAsFactor=F)
GEP<-t(GEP)
colnames(GEP)<-paste0(CELL, "_P", c(1:K))
if(CELL==cell_list[1]){
summary<-GEP
} else {
summary<-merge(summary, GEP, by=0, all=T)
row.names(summary)<-summary[,1]
summary<-summary[,-1]
}
}

##calculate jaccard index

tmp<-combn(ncol(summary), 2)
jaccard <- function(a, b) {
    intersection = length(intersect(a, b))
    union = length(a) + length(b) - intersection
    return (intersection/union)
}

jaccard_df<-data.frame()
for(i in 1:ncol(summary)){
for(j in 1:ncol(summary)){

A<-row.names(summary)[!is.na(summary[,i])][ order(summary[,i][!is.na(summary[,i])], decreasing=T)][1:300]
B<-row.names(summary)[!is.na(summary[,j])][ order(summary[,j][!is.na(summary[,j])], decreasing=T)][1:300]
JI=jaccard(A,B)
hoge<-data.frame(A=colnames(summary)[i], B=colnames(summary)[j], JI=JI)
jaccard_df<-rbind(jaccard_df, hoge)
}}

write.table(jaccard_df, "multiCell_model/Program_genes_Jaccard.txt", row.names=F, sep="\t", quote=F)
