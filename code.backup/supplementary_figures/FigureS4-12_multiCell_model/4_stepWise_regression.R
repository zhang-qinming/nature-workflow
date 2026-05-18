args   <- commandArgs(trailingOnly = T)  
options("scipen"=10) ##disable 100000 to 1e6
FILE<-as.character( args[1] ) ##e.g., Backman_2021_133.per_gene_estimates.tsv

K<-60 ##60

##data loading
cell_list=dir("data/Perturbseq/cNMF_regulation/")
K=60

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


##remove genes with na ..1942 genes left
summary<-na.omit(summary)
##
summary<-data.frame(gene=row.names(summary), summary)

corresp<-read.table("data/gencode_v41_gname_gid_ALL_sorted_onlyID", header=F, stringsAsFactor=F)
corresp<-corresp[!duplicated(corresp[,1]),]
corresp<-corresp[!duplicated(corresp[,2]),]
row.names(corresp)<-corresp[,1]

LOF<-read.table(paste0("data/LoF/GeneBayes_posterior/", FILE), sep="\t", quote="", header=T, stringsAsFactor=F)

LOF$post_mean[is.element(LOF$post_mean, "Inf")]<-max(LOF$post_mean[!is.infinite(LOF$post_mean)])
LOF$post_mean[is.element(LOF$post_mean, "-Inf")]<-min(LOF$post_mean[!is.infinite(LOF$post_mean)])

LOF<-data.frame(LOF, gene=corresp[as.character(LOF$ensg),2])

shet<-read.table("data/shet_10bins.txt", header=T, stringsAsFactor=F)
shet<-shet[is.element(shet$ensg, LOF$ensg),]


df<-merge(LOF, summary, by="gene")
df<-merge(df, shet, by="ensg")
df<-df[,c("post_mean", colnames(df)[grep("_P", colnames(df))], "shet")]

for(i in 2:ncol(df)){
df[,i][is.infinite(df[,i])] <-max(df[,i][!is.infinite(df[,i])])
}

##regress shet beforehand
model <- lm(post_mean ~ shet, data = df)
resids <- model$residuals
df$post_mean<-resids
df<-df[,!is.element(colnames(df), "shet")]


library(leaps)
b <- leaps::regsubsets(post_mean ~ ., data=df, nbest=1, nvmax=10, method="forward", really.big=T)

fit_all_sum = summary(b)
fit_all_sum_w=fit_all_sum[[1]]

pdf(paste0("multiCell_model/", FILE, "_BICetc.pdf"), width=10, height=10)
par(mfrow = c(2, 2))
plot(fit_all_sum$rss, xlab = "Number of Variables", ylab = "RSS", type = "b")

plot(fit_all_sum$adjr2, xlab = "Number of Variables", ylab = "Adjusted RSq", type = "b")
best_adj_r2 = which.max(fit_all_sum$adjr2)
points(best_adj_r2, fit_all_sum$adjr2[best_adj_r2],
       col = "red",cex = 2, pch = 20)

plot(fit_all_sum$cp, xlab = "Number of Variables", ylab = "Cp", type = 'b')
best_cp = which.min(fit_all_sum$cp)
points(best_cp, fit_all_sum$cp[best_cp], 
       col = "red", cex = 2, pch = 20)

plot(fit_all_sum$bic, xlab = "Number of Variables", ylab = "BIC", type = 'b')
best_bic = which.min(fit_all_sum$bic)
points(best_bic, fit_all_sum$bic[best_bic], 
       col = "red", cex = 2, pch = 20)

dev.off()

res<-data.frame(RSS=fit_all_sum$rss, adjr2=fit_all_sum$adjr2, Cp=fit_all_sum$cp, BIC=fit_all_sum$bic)

Programs=data.frame(Nprogram=c(1:5), Programs=c(paste(colnames(fit_all_sum_w)[fit_all_sum_w[1,]], collapse=","), paste(colnames(fit_all_sum_w)[fit_all_sum_w[2,]], collapse=","),
	paste(colnames(fit_all_sum_w)[fit_all_sum_w[3,]], collapse=","), paste(colnames(fit_all_sum_w)[fit_all_sum_w[4,]], collapse=","), paste(colnames(fit_all_sum_w)[fit_all_sum_w[5,]], collapse=",")))

write.table(res, paste0("multiCell_model/", FILE, "_res.txt"), row.names=F, sep="\t", quote=F)
write.table(Programs, paste0("multiCell_model/", FILE, "_selectedPrograms.txt"), row.names=F, sep="\t", quote=F)
