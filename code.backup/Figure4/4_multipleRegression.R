args   <- commandArgs(trailingOnly = T)  
options("scipen"=10)

K<-as.numeric( args[1] ) ##60

##data loading

GEP_reg_beta<-data.frame()

for(i in 1:K){
tmp<-read.table(paste0("data/Perturbseq/cNMF_regulation/K562GW/K",  K, "_program", i, "_perturb_effects.txt"), header=T)
tmp1<-tmp[,c(1,2)]
colnames(tmp1)<-c("GENE", paste0("P", i))

if(i==1){
GEP_reg_beta<-tmp1
} else {
GEP_reg_beta<-merge(GEP_reg_beta, tmp1, by="GENE")
}
}

##
corresp<-read.table("data/gencode_v41_gname_gid_ALL_sorted_onlyID", header=F, stringsAsFactor=F)
corresp<-corresp[!duplicated(corresp[,1]),]
row.names(corresp)<-corresp[,1]

##MCH
MCH<-read.table("data/LoF/GeneBayes_posterior/Backman_2021_86.per_gene_estimates.tsv", sep="\t", quote="", header=T, stringsAsFactor=F)
MCH<-data.frame(GENE=corresp[MCH$ensg, 2], MCH)

##RDW
RDW<-read.table("data/LoF/GeneBayes_posterior/Backman_2021_88.per_gene_estimates.tsv", sep="\t", quote="", header=T, stringsAsFactor=F)
RDW<-data.frame(GENE=corresp[RDW$ensg, 2], RDW)

##Immature ret fra
IRF<-read.table("data/LoF/GeneBayes_posterior/Backman_2021_106.per_gene_estimates.tsv", sep="\t", quote="", header=T, stringsAsFactor=F)
IRF<-data.frame(GENE=corresp[IRF$ensg, 2], IRF)


shet<-read.table("data/shet_10bins.txt", header=T, stringsAsFactor=F)
shet<-shet[is.element(shet$ensg, MCH$ensg),]
shet<-shet[,c("GENE", "shet")]

##1. MCH, P25 + P16 + P4
df<-merge(GEP_reg_beta, MCH, by="GENE")
df<-merge(df, shet, by="GENE")
df$post_mean<-scale(df$post_mean)
df$shet<-scale(df$shet)
df$P25<-scale(df$P25)
df$P16<-scale(df$P16)
df$P4<-scale(df$P4)

fit<-lm(post_mean ~ shet  + P25 + P16 + P4, data=df)

res<-as.data.frame(summary(fit)$coefficients)
library(ggplot2)
res<-res[-1,]
res<-res[order(abs(res$Estimate), decreasing=T),]
res<-data.frame(covariate=row.names(res), res)
res$covariate<-factor(as.character(res$covariate), levels=rev(as.character(res$covariate)))
colnames(res)[3]<-"SE"
res1=res[-1,]
res1$covariate<-gsub("P4", "S phase", res1$covariate)
res1$covariate<-gsub("P16", "Autophagy", res1$covariate)
res1$covariate<-gsub("P25", "G2M phase", res1$covariate)

g<-ggplot(res1, aes(x=covariate, y=Estimate))
g<-g +  theme_classic(base_size = 40, base_family = "Helvetica")
g<-g + geom_point(size=7)
g<-g + geom_errorbar(aes(ymax=Estimate + 1.96*SE, ymin=Estimate - 1.96*SE), width=0)
g<-g + geom_hline(yintercept=0, linetype="dashed")
g<-g + ylab("Standardized effect size") + xlab("")
g<-g + coord_flip()
ggsave(plot=g, paste0("Fig5F.pdf"), width=12, height=6 )
