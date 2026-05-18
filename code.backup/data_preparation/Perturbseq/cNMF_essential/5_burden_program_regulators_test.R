args   <- commandArgs(trailingOnly = T)  
options("scipen"=10) ##disable 100000 to 1e6
FILE<-as.character( args[1] ) ##e.g., Backman_2021_86.per_gene_estimates.tsv
K<-as.numeric( args[2] ) ##60

##data loading
GEP<-read.table(paste0("data/Perturbseq/cNMF/", CELL, "/test1/test1.gene_spectra_score.k_", K, ".dt_0_5.txt"), header=T, stringsAsFactor=F)
GEP<-t(GEP)
colnames(GEP)<-paste0("P", c(1:K))
corresp<-read.table("data/gencode_v41_gname_gid_ALL_sorted_onlyID", header=F, stringsAsFactor=F)
corresp<-corresp[!duplicated(corresp[,1]),]
row.names(corresp)<-corresp[,1]
GEP<-data.frame(GENE=corresp[row.names(GEP), 2], GEP)

GEP_reg_beta<-data.frame()

for(i in 1:K){
tmp<-read.table(paste0("data/Perturbseq/cNMF_regulation/", CELL, "/K", K, "_program", i, "_perturb_effects.txt"), header=T)
tmp1<-tmp[,c(1,2)]

colnames(tmp1)<-c("GENE", paste0("P", i))

if(i==1){
GEP_reg_beta<-tmp1
} else {
GEP_reg_beta<-merge(GEP_reg_beta, tmp1, by="GENE")
}
}

LOF<-read.table(paste0("data/LoF/GeneBayes_posterior/", FILE), sep="\t", quote="", header=T, stringsAsFactor=F)

LOF$post_mean[is.element(LOF$post_mean, "Inf")]<-max(LOF$post_mean[!is.infinite(LOF$post_mean)])
LOF$post_mean[is.element(LOF$post_mean, "-Inf")]<-min(LOF$post_mean[!is.infinite(LOF$post_mean)])

LOF<-data.frame(LOF, gene=corresp[as.character(LOF$ensg),2])

GEP2<-data.frame(ensg=row.names(GEP), GEP)
df<-merge(GEP2, LOF, by="ensg")


shet<-read.table("data/shet_10bins.txt", header=T, stringsAsFactor=F)
shet<-shet[is.element(shet$ensg, LOF$ensg),]

summary_reg<-data.frame()
summary_pro<-data.frame()

for(Program in c(1:K)){

tmp<-GEP_reg_beta[,c("GENE", paste0("P", Program))]
colnames(tmp)<-c("GENE", "perturb_beta")
tmp<-data.frame(tmp, ensg=corresp[ as.character(tmp$GENE), 1])
df<-merge(tmp, LOF, by="ensg")

P_pearson_all<-cor.test(df$perturb_beta, df$post_mean, method="pearson")$p.value
R_pearson_all<-cor.test(df$perturb_beta, df$post_mean, method="pearson")$estimate
R_pearson_all_CIlower<-cor.test(df$perturb_beta, df$post_mean, method="pearson")$conf.int[1]
R_pearson_all_CIupper<-cor.test(df$perturb_beta, df$post_mean, method="pearson")$conf.int[2]

##with shet regression
df<-merge(df, shet, by="ensg")
df$post_mean<-scale(df$post_mean)
df$perturb_beta<-scale(df$perturb_beta)

fit<- lm(post_mean~perturb_beta + shet, data=df)
P_withShet<-summary(fit)$coefficients[2,4]
beta_withShet<-summary(fit)$coefficients[2,1]
betaSE_withShet<-summary(fit)$coefficients[2,2]


hoge<-data.frame(FILE=FILE, Program=Program, 
    P_pearson_all=P_pearson_all, R_pearson_all=R_pearson_all,
    R_pearson_all_CIlower=R_pearson_all_CIlower, R_pearson_all_CIupper=R_pearson_all_CIupper,
    P_withShet=P_withShet, beta_withShet=beta_withShet, betaSE_withShet=betaSE_withShet)
summary_reg<-rbind(summary_reg, hoge)

###
##from here, program enrichment

##1. mean(gamma) for top 100 program genes

tmp<-GEP[,c("GENE",paste0("P", Program))]
tmp<-tmp[order(tmp[,2], decreasing=T),]
colnames(tmp)[2]<-"GEPscore"
tmp<-data.frame(tmp, ensg=row.names(tmp))
df<-merge(tmp, LOF, by="ensg")
df2<-df[order(df$GEPscore, decreasing=T),][1:100,]

MEANgamma_top100=mean(df2[,"post_mean"])

 
shet_tmp<-shet[is.element(shet$ensg, df2$ensg),]
shet_tmp_all<-shet[is.element(shet$ensg, df$ensg),]

A<-table(shet_tmp$shet_BIN)

random<-c()
for(I in 1:100000){
set.seed(I)
genes<-c()
for(p in 1:length(A)){
genes<-c(genes, sample(shet_tmp_all$ensg[is.element(shet_tmp_all$shet_BIN, names(A)[p])], A[p]))
}
random<-c(random, mean(df[is.element(df$ensg, genes),"post_mean"], na.omit=T))
}

random_mean=mean(random)


P1=(rank(c(MEANgamma_top100, random))[1]/length(random)) * 2
P2=(( length(random) + 2 - (rank( c(MEANgamma_top100, random))[1]) ) /length(random)) * 2
program_P=min(P1,P2)


hoge<-data.frame(FILE=FILE, Program=Program, 
MEANgamma_top100=MEANgamma_top100,
shet_adjusted_random_mean=random_mean,
MEANgamma_top100_shet_adjusted_P=program_P
)

summary_pro<-rbind(summary_pro, hoge)

}

write.table(summary_reg, paste0("data/Perturbseq/trait_association/", CELL, "/ProgramLevel/regulators_enrichment_K", K, "_", FILE), row.names=F, sep="\t", quote=F)
write.table(summary_pro, paste0("data/Perturbseq/trait_association/", CELL, "/ProgramLevel/programs_enrichment_K", K, "_", FILE), row.names=F, sep="\t", quote=F)
