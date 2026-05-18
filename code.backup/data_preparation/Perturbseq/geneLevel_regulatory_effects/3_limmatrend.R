args   <- commandArgs(trailingOnly = T)   # enable reading factors from linux command line

options("scipen"=10) ##disable 100000 to 1e6
dir.create("gwps_limma", showWarnings=F)

AA<-as.numeric( args[1] )

library(dplyr)
library(kSamples)

library(data.table)

dfA<-fread(paste0("limmatrend/tmp_limma_input/", AA, ".input"), header=T, data.table=F)
dfB<-fread(paste0("limmatrend/tmp_limma_input/control.input"), header=T, data.table=F)
df_all<-rbind(dfA, dfB)

N=ncol(df_all)-6

GENES<-unique(df_all$gene)
GENES<-GENES[!is.element(GENES, "non-targeting")]

library(limma)
library(edgeR)
for(GENE in GENES){
if(!file.exists(paste0("gwps_limma/", GENE,  "_KO.txt"))){
df<-df_all[is.element(df_all$gene, c(GENE, "non-targeting")),]

df$gem_group <-factor(as.character(df$gem_group))

exp<-t(df[,2:N])
colnames(exp)<-df[,1]
COND<-cbind(model.matrix(~GROUP + gem_group, df),df[,c("n_genes", "mitopercent")] )

dge <- DGEList(exp, group = COND$GROUPperturb)
dge <- calcNormFactors(dge)
y <- new("EList")
y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 2)

fit <- lmFit(y, design=COND)
fit <- eBayes(fit, trend = TRUE, robust = TRUE)

result <- as.data.frame(topTable(fit, coef="GROUPperturb", sort.by = "P", n = Inf))

write.table(data.frame(GENE=row.names(result), result), paste0("gwps_limma/", GENE,  "_KO.txt"), row.names=F, sep="\t", quote=F)
}}
