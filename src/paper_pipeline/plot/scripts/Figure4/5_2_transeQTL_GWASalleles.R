args   <- commandArgs(trailingOnly = T)  
options("scipen"=10) ##disable 100000 to 1e6
P<-as.character( args[1] ) ##1-60
topN <- if (length(args) >= 2) as.numeric(args[2]) else 100
gwas_path <- if (length(args) >= 3) as.character(args[3]) else "data/GWAS/30050_irnt.tsv.gz"
trait_label <- if (length(args) >= 4) as.character(args[4]) else "MCH"
dir.create("program_trans/result", showWarnings=F)

library(data.table)
MCH_GWAS<-fread(gwas_path, header=T, data.table=F)


MCH_GWAS_ctr<-MCH_GWAS[MCH_GWAS$pval>0.05,]
MCH_GWAS_ctr<-na.omit(MCH_GWAS_ctr)
MCH_GWAS_ctr<-data.frame(MCH_GWAS_ctr, CHR=sapply(strsplit(as.character(MCH_GWAS_ctr$variant), ":"), function(x){x[1]}), POS=sapply(strsplit(as.character(MCH_GWAS_ctr$variant), ":"), function(x){x[2]}), REF=sapply(strsplit(as.character(MCH_GWAS_ctr$variant), ":"), function(x){x[3]}), ALT=sapply(strsplit(as.character(MCH_GWAS_ctr$variant), ":"), function(x){x[4]}))
MCH_GWAS_ctr<-data.frame(MCH_GWAS_ctr, ID=paste0(MCH_GWAS_ctr$CHR, ":", MCH_GWAS_ctr$POS))


MCH_GWAS<-MCH_GWAS[MCH_GWAS$pval<5e-8,]
MCH_GWAS<-na.omit(MCH_GWAS)
MCH_GWAS<-data.frame(MCH_GWAS, CHR=sapply(strsplit(as.character(MCH_GWAS$variant), ":"), function(x){x[1]}), POS=sapply(strsplit(as.character(MCH_GWAS$variant), ":"), function(x){x[2]}), REF=sapply(strsplit(as.character(MCH_GWAS$variant), ":"), function(x){x[3]}), ALT=sapply(strsplit(as.character(MCH_GWAS$variant), ":"), function(x){x[4]}))
MCH_GWAS<-data.frame(MCH_GWAS, ID=paste0(MCH_GWAS$CHR, ":", MCH_GWAS$POS))

##

res_summary<-data.frame()

data<-read.table(paste0("program_trans/P", P, "_top", topN, "genes_trans.txt"), header=T, stringsAsFactor=F)
data<-data.frame(data, ID=paste0(data$SNPChr, ":", data$SNPPos))

data_f<-merge(MCH_GWAS, data, by="ID")
summary<-data.frame()
for(i in 1:nrow(data_f)){
if(data_f[i,"ALT"] == data_f[i,"AssessedAllele"]){
Z_risk=sign(data_f[i,"beta"]) * data_f[i,"Zscore"]
hoge<-data.frame(ID=data_f$ID[i], SNP=data_f$SNP[i], GENE=data_f$GeneSymbol[i], ensg=data_f$Gene[i], Z_risk=Z_risk)
summary<-rbind(summary, hoge)
} else if (data_f[i,"REF"] == data_f[i,"AssessedAllele"]){
Z_risk= -sign(data_f[i,"beta"]) * data_f[i,"Zscore"]
hoge<-data.frame(ID=data_f$ID[i], SNP=data_f$SNP[i], GENE=data_f$GeneSymbol[i], ensg=data_f$Gene[i], Z_risk=Z_risk)
summary<-rbind(summary, hoge)
}
}


data_ff<-merge(MCH_GWAS_ctr, data, by="ID")
summary_ctr<-data.frame()
for(i in 1:nrow(data_ff)){
if(data_ff[i,"ALT"] == data_ff[i,"AssessedAllele"]){
Z_risk=sign(data_ff[i,"beta"]) * data_ff[i,"Zscore"]
hoge<-data.frame(ID=data_ff$ID[i], SNP=data_ff$SNP[i], GENE=data_ff$GeneSymbol[i], ensg=data_ff$Gene[i], Z_risk=Z_risk)
summary_ctr<-rbind(summary_ctr, hoge)
} else if (data_ff[i,"REF"] == data_ff[i,"AssessedAllele"]){
Z_risk= -sign(data_ff[i,"beta"]) * data_ff[i,"Zscore"]
hoge<-data.frame(ID=data_ff$ID[i], SNP=data_ff$SNP[i], GENE=data_ff$GeneSymbol[i], ensg=data_ff$Gene[i], Z_risk=Z_risk)
summary_ctr<-rbind(summary_ctr, hoge)
}
}


##take mean
library(plyr)
summary2 <- ddply(.data=summary, 
                 .(SNP), 
                 summarize, 
                Z_risk=mean(Z_risk))

summary_ctr2<-ddply(.data=summary_ctr, 
                 .(SNP), 
                 summarize, 
                Z_risk=mean(Z_risk))


##T test
P_T=t.test(summary_ctr2$Z_risk, summary2$Z_risk)$p.value

mean_ctr=mean(summary_ctr2$Z_risk)
mean_GWAS=mean(summary2$Z_risk)

hoge<-data.frame(Program=P, topN=topN, Ttest_P=P_T, mean_Z_GWAS=mean_GWAS, mean_Z_ctrl=mean_ctr)

res_summary<-rbind(res_summary, hoge)



write.table(
  res_summary,
  paste0("program_trans/result/P", P, "_transeQTL_", trait_label, "_directionalTest.txt"),
  row.names=F,
  sep="\t",
  quote=F
)
