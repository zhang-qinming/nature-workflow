
args <- commandArgs(trailingOnly = TRUE)
K <- if (length(args) >= 1) as.numeric(args[1]) else 60
mch_posterior_path <- if (length(args) >= 2) as.character(args[2]) else "data/LoF/GeneBayes_posterior/Backman_2021_86.per_gene_estimates.tsv"
rdw_posterior_path <- if (length(args) >= 3) as.character(args[3]) else "data/LoF/GeneBayes_posterior/Backman_2021_88.per_gene_estimates.tsv"
irf_posterior_path <- if (length(args) >= 4) as.character(args[4]) else "data/LoF/GeneBayes_posterior/Backman_2021_106.per_gene_estimates.tsv"

##data loading
GEP<-read.table(paste0("data/Perturbseq/cNMF/K562GW/cNMF_all.gene_spectra_score.k_", K, ".dt_0_5.txt"), header=T, stringsAsFactor=F)
GEP<-t(GEP)
colnames(GEP)<-paste0("P", c(1:K))
corresp<-read.table("data/gencode_v41_gname_gid_ALL_sorted_onlyID", header=F, stringsAsFactor=F)
corresp<-corresp[!duplicated(corresp[,1]),]
row.names(corresp)<-corresp[,1]
GEP<-data.frame(GENE=corresp[row.names(GEP), 2], GEP)

GEP_reg_beta<-data.frame()
GEP_reg_P<-data.frame()

for(i in 1:K){
tmp<-read.table(paste0("data/Perturbseq/cNMF_regulation/K562GW/K",  K, "_program", i, "_perturb_effects.txt"), header=T)

tmp1<-tmp[,c(1,2)]
tmp2<-tmp[,c(1,3)]

colnames(tmp1)<-c("GENE", paste0("P", i))
colnames(tmp2)<-c("GENE", paste0("P", i))

if(i==1){
GEP_reg_beta<-tmp1
GEP_reg_P<-tmp2
} else {
GEP_reg_beta<-merge(GEP_reg_beta, tmp1, by="GENE")
GEP_reg_P<-merge(GEP_reg_P, tmp2, by="GENE")
}
}

##MCH
MCH<-read.table(mch_posterior_path, sep="\t", quote="", header=T, stringsAsFactor=F)
MCH<-data.frame(gene=corresp[MCH$ensg, 2], MCH)

##RDW
RDW<-read.table(rdw_posterior_path, sep="\t", quote="", header=T, stringsAsFactor=F)
RDW<-data.frame(gene=corresp[RDW$ensg, 2], RDW)

##Immature ret fra
IRF<-read.table(irf_posterior_path, sep="\t", quote="", header=T, stringsAsFactor=F)
IRF<-data.frame(gene=corresp[IRF$ensg, 2], IRF)


##
shet<-read.table("data/shet_10bins.txt", header=T, stringsAsFactor=F)
shet<-shet[is.element(shet$ensg, MCH$ensg),]
shet<-shet[,c("ensg", "shet")]

df<-merge(MCH, shet, by="ensg")
model<-lm(post_mean~shet, data=df)
resid<-model$residuals
MCH<-data.frame(resid=resid, df)


df<-merge(RDW, shet, by="ensg")
model<-lm(post_mean~shet, data=df)
resid<-model$residuals
RDW<-data.frame(resid=resid, df)


df<-merge(IRF, shet, by="ensg")
model<-lm(post_mean~shet, data=df)
resid<-model$residuals
IRF<-data.frame(resid=resid, df)

summary<-data.frame()

for(SEED in c(1:500)){

set.seed(SEED)
MCH_train<-MCH[sample(nrow(MCH), round(0.8*nrow(MCH))),]
MCH_test<-MCH[!is.element(MCH$ensg, MCH_train$ensg),]

df<-merge(MCH_train, GEP_reg_beta, by="GENE")

fit<-lm(resid ~ P40 + P25 + P16, data=df)

df2<-merge(MCH_test, GEP_reg_beta, by="GENE")
predict=predict(fit,df2) #Predictions on Testing data

RSQ = cor(df2$resid, predict)^2
adjRsq = 1 - ((1 - RSQ) * (length(predict) - 1) / (length(predict) - 3 - 1))

hoge<-data.frame(TRAIT="MCH", model="3program_selected", SEED=SEED, RSQ=RSQ, adjRsq=adjRsq)
summary<-rbind(summary, hoge)

df_f<-df[,c(2, sample(12:ncol(df), 3) )]
fit<-lm(resid ~ ., data=df_f)
predict=predict(fit,df2) #Predictions on Testing data
RSQ = cor(df2$resid, predict)^2
adjRsq = 1 - ((1 - RSQ) * (length(predict) - 1) / (length(predict) - 3 - 1))

hoge<-data.frame(TRAIT="MCH", model="3program_random", SEED=SEED, RSQ=RSQ, adjRsq=adjRsq)
summary<-rbind(summary, hoge)

for(P in seq_len(K)){
df_f<-df[,c("resid", paste0("P", P) )]
fit<-lm(resid ~ ., data=df_f)
predict=predict(fit,df2) #Predictions on Testing data
RSQ = cor(df2$resid, predict)^2
adjRsq = 1 - ((1 - RSQ) * (length(predict) - 1) / (length(predict) - 1 - 1))
hoge<-data.frame(TRAIT="MCH", model=paste0("P", P), SEED=SEED, RSQ=RSQ, adjRsq=adjRsq)
summary<-rbind(summary, hoge)

}

set.seed(SEED)
RDW_train<-RDW[sample(nrow(RDW), round(0.8*nrow(RDW))),]
RDW_test<-RDW[!is.element(RDW$ensg, RDW_train$ensg),]

df<-merge(RDW_train, GEP_reg_beta, by="GENE")

fit<-lm(resid ~ P16 + P43 + P25, data=df)

df2<-merge(RDW_test, GEP_reg_beta, by="GENE")
predict=predict(fit,df2) #Predictions on Testing data

RSQ = cor(df2$resid, predict)^2
adjRsq = 1 - ((1 - RSQ) * (length(predict) - 1) / (length(predict) - 3 - 1))

hoge<-data.frame(TRAIT="RDW", model="3program_selected", SEED=SEED, RSQ=RSQ, adjRsq=adjRsq)
summary<-rbind(summary, hoge)

df_f<-df[,c(2, sample(12:ncol(df), 3) )]
fit<-lm(resid ~ ., data=df_f)
predict=predict(fit,df2) #Predictions on Testing data
RSQ = cor(df2$resid, predict)^2
adjRsq = 1 - ((1 - RSQ) * (length(predict) - 1) / (length(predict) - 3 - 1))

hoge<-data.frame(TRAIT="RDW", model="3program_random", SEED=SEED, RSQ=RSQ, adjRsq=adjRsq)
summary<-rbind(summary, hoge)

for(P in seq_len(K)){
df_f<-df[,c("resid", paste0("P", P) )]
fit<-lm(resid ~ ., data=df_f)
predict=predict(fit,df2) #Predictions on Testing data
RSQ = cor(df2$resid, predict)^2
adjRsq = 1 - ((1 - RSQ) * (length(predict) - 1) / (length(predict) - 1 - 1))

hoge<-data.frame(TRAIT="RDW", model=paste0("P", P), SEED=SEED, RSQ=RSQ, adjRsq=adjRsq)
summary<-rbind(summary, hoge)

}
}

##figure

res<-summary[is.element(summary$TRAIT, "MCH"),]
library(plyr)
tmp <- ddply(.data=res, 
                 .(model), 
                 summarize, 
               M=median(RSQ))

tmp2<-tmp[order(tmp$M, decreasing=T),]
selected<-tmp2$model[grep("P", tmp2$model)][1:5]

mat<-res[is.element(res$model, c("3program_selected", "3program_random", selected)),]

library("plotrix")

se_sum <- ddply(.data=mat, 
                 .(model), 
                 summarize,
                 median=median(RSQ),
               ymin=median(RSQ) - 2*std.error(RSQ), ymax=median(RSQ) + 2*std.error(RSQ) )


se_sum$model<-factor(as.character(se_sum$model), levels=c("3program_selected", "3program_random", selected))

library(ggplot2)

g<-ggplot(se_sum, aes(x=model, y=median))
g<-g +  theme_classic(base_size = 20, base_family = "Helvetica")
g<-g + geom_bar(stat="identity", position="dodge")
g<-g + geom_errorbar(data=se_sum, aes(ymin = ymin, ymax = ymax), position = "dodge", linewidth=1, width=0)
g<-g + ylab("r squared")
g<-g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot=g, "FigS7_MCH_model_rsquared.pdf", width=12, height=6)

##

res<-summary[is.element(summary$TRAIT, "RDW"),]
library(plyr)
tmp <- ddply(.data=res, 
                 .(model), 
                 summarize, 
               M=median(RSQ))

tmp2<-tmp[order(tmp$M, decreasing=T),]
selected<-tmp2$model[grep("P", tmp2$model)][1:5]

mat<-res[is.element(res$model, c("3program_selected", "3program_random", selected)),]

library("plotrix")

se_sum <- ddply(.data=mat, 
                 .(model), 
                 summarize,
                 median=median(RSQ),
               ymin=median(RSQ) - 2*std.error(RSQ), ymax=median(RSQ) + 2*std.error(RSQ) )

se_sum$model<-factor(as.character(se_sum$model), levels=c("3program_selected", "3program_random", selected))

library(ggplot2)

g<-ggplot(se_sum, aes(x=model, y=median))
g<-g +  theme_classic(base_size = 20, base_family = "Helvetica")
g<-g + geom_bar(stat="identity", position="dodge")
g<-g + geom_errorbar(data=se_sum, aes(ymin = ymin, ymax = ymax), position = "dodge", linewidth=1, width=0)
g<-g + ylab("r squared")
g<-g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot=g, "FigS7_RDW_model_rsquared.pdf", width=12, height=6)
