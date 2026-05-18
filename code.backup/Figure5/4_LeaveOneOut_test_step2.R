
args   <- commandArgs(trailingOnly = T)   # enable reading factors from linux command line
options("scipen"=10) ##disable 100000 to 1e6
genes<-as.character( args[1] ) ##file with -100 gene names. split genes into chunks before running this code.

Program_N=as.numeric( args[2] ) ##5
Regulator_N=as.numeric( args[3] ) ##3
TRAIT=as.character( args[4] ) ##"MCH"
Program_top_def=as.numeric( args[5] ) ##200

K=60

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
MCH<-read.table("data/LoF/GeneBayes_posterior/Backman_2021_86.per_gene_estimates.tsv", sep="\t", quote="", header=T, stringsAsFactor=F)
MCH<-data.frame(gene=corresp[MCH$ensg, 2], MCH)

##RDW
RDW<-read.table("data/LoF/GeneBayes_posterior/Backman_2021_88.per_gene_estimates.tsv", sep="\t", quote="", header=T, stringsAsFactor=F)
RDW<-data.frame(gene=corresp[RDW$ensg, 2], RDW)

##Immature ret fra
IRF<-read.table("data/LoF/GeneBayes_posterior/Backman_2021_106.per_gene_estimates.tsv", sep="\t", quote="", header=T, stringsAsFactor=F)
IRF<-data.frame(gene=corresp[IRF$ensg, 2], IRF)


##
shet<-read.table("data/shet_10bins.txt", header=T, stringsAsFactor=F)
shet<-shet[is.element(shet$ensg, MCH$ensg),]

##
shet<-read.table("data/shet_10bins.txt", header=T, stringsAsFactor=F)
shet<-shet[is.element(shet$ensg, MCH$ensg),]


gene_list<-read.table(paste0( "loo_topgene_prediction/gene_list/", genes ), header=F, stringsAsFactor=F)[,1]

summary_all<-data.frame()
for(loo_gene in gene_list){

##
if(TRAIT=="MCH"){
LOF=MCH
}
if(TRAIT=="RDW"){
LOF=RDW
}
if(TRAIT=="IRF"){
LOF=IRF
}

##leave one out
gamma_loo=LOF[is.element(LOF$GENE, loo_gene),"post_mean"]

LOF<-LOF[!is.element(LOF$GENE, loo_gene),]

##

summary<-data.frame()

##step1. estimate program burden effects to find the defined number of top enriched programs

Program_P_sum<-data.frame()
for(Program in 1:K){

tmp<-GEP[,c("GENE",paste0("P", Program))]
tmp<-tmp[order(tmp[,2], decreasing=T),]
colnames(tmp)[2]<-"GEPscore"
tmp<-data.frame(tmp, ensg=row.names(tmp))
df<-merge(tmp, LOF, by="ensg")
df2<-df[order(df$GEPscore, decreasing=T),][1:Program_top_def,]

MEANgamma_top100=mean(df2[,"post_mean"])

 
##make random 10000 set of genes, matched to shet bin

shet_tmp<-shet[is.element(shet$ensg, df2$ensg),]
shet_tmp_all<-shet[is.element(shet$ensg, df$ensg),]

A<-table(shet_tmp$shet_BIN)

random<-c()
for(I in 1:10000){
set.seed(I)
genes2<-c()
for(p in 1:length(A)){
genes2<-c(genes2, sample(shet_tmp_all$ensg[is.element(shet_tmp_all$shet_BIN, names(A)[p])], A[p]))
}
random<-c(random, mean(df[is.element(df$ensg, genes2),"post_mean"], na.omit=T))
}

random_mean=mean(random)


P1=(rank(c(MEANgamma_top100, random))[1]/length(random)) * 2
P2=(( length(random) + 2 - (rank( c(MEANgamma_top100, random))[1]) ) /length(random)) * 2
program_P=min(P1,P2)

hoge<-data.frame(Program=Program, P=program_P, meanG=MEANgamma_top100 - random_mean )
Program_P_sum=rbind(Program_P_sum, hoge)
}
Program_P_sum<-Program_P_sum[order(abs(Program_P_sum$meanG), decreasing=T),]
Program_P_sum<-Program_P_sum[order(Program_P_sum$P),]

Program_selected=paste0("P", Program_P_sum$Program[1:Program_N])


##step2. test regulator-burden correlations to find the defined number of top enriched programs

df<-merge(LOF, GEP_reg_beta, by="GENE")
df<-merge(df, shet, by="ensg")
df<-df[,c("post_mean", paste0("P", c(1:60)), "shet")]
for(i in 2:ncol(df)){
df[,i][is.infinite(df[,i])] <-max(df[,i][!is.infinite(df[,i])])
}

library(leaps)
b <- leaps::regsubsets(post_mean ~ ., data=df, nbest=1, nvmax=Regulator_N+1,really.big=T)
fit_all_sum = summary(b)
fit_all_sum_w=fit_all_sum[[1]]
regulator_selected=colnames(fit_all_sum_w)[fit_all_sum_w[nrow(fit_all_sum_w),]]
regulator_selected=regulator_selected[!is.element(regulator_selected, c("(Intercept)", "shet"))]

if(length(regulator_selected) !=Regulator_N ){
b <- leaps::regsubsets(post_mean ~ ., data=df, nbest=1, nvmax=Regulator_N,really.big=T)
fit_all_sum = summary(b)
fit_all_sum_w=fit_all_sum[[1]]
regulator_selected=colnames(fit_all_sum_w)[fit_all_sum_w[nrow(fit_all_sum_w),]]
regulator_selected=regulator_selected[!is.element(regulator_selected, c("(Intercept)", "shet"))]
}

##step3. predict the direction of associations of genes in the model

df<-merge(LOF, GEP_reg_beta, by="GENE")
df<-merge(df, shet, by="ensg")
df<-df[,c("GENE.x", "post_mean", regulator_selected, "shet")]
colnames(df)[1]<-"GENE"
for(i in 2:ncol(df)){
df[,i][is.infinite(df[,i])] <-max(df[,i][!is.infinite(df[,i])])
df[,i]<-scale(df[,i])
}


fit1<-lm(post_mean~., data=df[,-1])
res1<-summary(fit1)$coefficients

df2<-GEP_reg_P
row.names(df2)<-df2$GENE

df3<-GEP_reg_beta
sum_effect<-data.frame(GENE=df3$GENE, effect=0, regeffect=0)
for(p in regulator_selected){
tmp<-res1[p, "Estimate"] *df3[,p]
tmp[p.adjust(df2[,p], method="BH")>0.05]<-0
sum_effect$effect=sum_effect$effect + tmp
}

sum_effect$regeffect<-sign(sum_effect$effect)
sum_effect2<-data.frame(GENE=GEP$GENE, proeffect=0)

for(p in rev(Program_selected)){  
tmp<-GEP[order(GEP[,p], decreasing=T), "GENE"][1:Program_top_def]

tmp2<-sign( Program_P_sum$meanG [Program_P_sum$Program==gsub("P", "", p)])
sum_effect2$proeffect[is.element(sum_effect2$GENE, tmp)]<-tmp2
}

sum_effect3<-merge(sum_effect, sum_effect2, by="GENE", all=TRUE)

sum_effect3[is.na(sum_effect3)]<-0

##
hoge<-data.frame(gamma_loo=gamma_loo, sum_effect3[is.element(sum_effect3$GENE, loo_gene),])

if(hoge$proeffect !=0 ){
	hoge$effect <-sign(hoge$proeffect)
} else {

	hoge$effect <-sign(hoge$regeffect)
}

summary_all<-rbind(summary_all, hoge)
}
write.table(summary_all, paste0("loo_topgene_prediction/", TRAIT, "_", genes, "_sign_predicted.txt"), row.names=F, sep="\t", quote=F, append=F)


##after finished, compare large gamma vs low gamma genes for the accuracy of prediction.
