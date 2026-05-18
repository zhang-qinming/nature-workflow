cell_list=dir("data/Perturbseq/cNMF_regulation")
K=60
dir.create("multiCell_model", showWarnings=F)

summary<-data.frame()
for(CELL in cell_list){
GEP_reg_beta<-data.frame()
for(i in 1:K){
tmp<-read.table(paste0("data/Perturbseq/cNMF_regulation/", CELL, "/K60_program", i, "_perturb_effects.txt"), header=T)
tmp1<-tmp[,c(1,2)]
colnames(tmp1)<-c("GENE", paste0("P", i))
if(i==1){
GEP_reg_beta<-tmp1
} else {
GEP_reg_beta<-merge(GEP_reg_beta, tmp1, by="GENE")
}
}
colnames(GEP_reg_beta)<-paste0(CELL, "_", colnames(GEP_reg_beta))
row.names(GEP_reg_beta)<-GEP_reg_beta[,1]
GEP_reg_beta<-GEP_reg_beta[,-1]
if(CELL==cell_list[1]){
summary<-GEP_reg_beta
} else {
summary<-merge(summary, GEP_reg_beta, by=0, all=T)
row.names(summary)<-summary[,1]
summary<-summary[,-1]
}}


##calculate correlation of regulatory effects

reg_df<-data.frame()
for(i in 1:ncol(summary)){
for(j in 1:ncol(summary)){
if(i !=j ){
df<-na.omit(summary[,c(i,j)])
P=cor.test(df[,1], df[,2])$p.value
R=cor.test(df[,1], df[,2])$estimate

hoge<-data.frame(A=colnames(summary)[i], B=colnames(summary)[j], P=P, R=R)
reg_df<-rbind(reg_df, hoge)
}}}

write.table(reg_df, "multiCell_model/Regulators_correlation.txt", row.names=F, sep="\t", quote=F)
