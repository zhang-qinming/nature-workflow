file_list=dir("gwps_limma")
file_list=file_list[grep("_KO.txt", file_list)]

for(FILE in file_list){
data<-read.table(paste0("gwps_limma/", FILE), header=T, stringsAsFactor=F)
GENE=gsub("_KO.txt", "", FILE)
if(FILE == file_list[1]){
beta_mat=data[,c("GENE", "logFC")]
colnames(beta_mat)[2]<-GENE
p_mat<-data[,c("GENE", "P.Value")]
colnames(p_mat)[2]<-GENE
} else {
beta_mat<-merge(beta_mat, data[,c("GENE", "logFC")], by="GENE")
colnames(beta_mat)[ncol(beta_mat)]<-GENE
p_mat<-merge(p_mat, data[,c("GENE", "P.Value")], by="GENE")
colnames(p_mat)[ncol(p_mat)]<-GENE
}
}

dir.create("data/Perturbseq/geneLevel/K562GW/", showWarnings=F)

write.table(beta_mat, "data/Perturbseq/geneLevel/K562GW/limma_logFC_sum.txt", row.names=F, sep="\t", quote=F)
write.table(p_mat, "data/Perturbseq/geneLevel/K562GW/limma_P_sum.txt", row.names=F, sep="\t", quote=F)
