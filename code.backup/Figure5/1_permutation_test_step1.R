
args   <- commandArgs(trailingOnly = T)   # enable reading factors from linux command line
options("scipen"=10) ##disable 100000 to 1e6
dir.create("permutation_test", showWarnings=F)
dir.create(paste0("permutation_test/P", Program_top_def), showWarnings=F)

Program_N=as.numeric( args[1] ) ##1-6
Regulator_N=as.numeric( args[2] ) ##1-6
Program_top_def=as.numeric( args[3] ) ##100, 200,300
LOF_thresh=as.numeric( args[4] ) ##0.03,0.05,0.1,0.12,0.14

TRAIT="MCH"
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
if(TRAIT=="MCH"){
LOF=MCH
}
if(TRAIT=="RDW"){
LOF=RDW
}
if(TRAIT=="IRF"){
LOF=IRF
}

#####
##analyze the observed data
#####

summary<-data.frame()

##step1. test program burden effects to find the defined number of top enriched programs

Program_P_sum<-data.frame()
for(Program in 1:K){

tmp<-GEP[,c("GENE",paste0("P", Program))]
tmp<-tmp[order(tmp[,2], decreasing=T),]
colnames(tmp)[2]<-"GEPscore"
tmp<-data.frame(tmp, ensg=row.names(tmp))
df<-merge(tmp, LOF, by="ensg")
df2<-df[order(df$GEPscore, decreasing=T),][1:Program_top_def,]

MEANgamma_top100=mean(df2[,"post_mean"] )

##make random 10000 set of genes, matched to shet bin

shet_tmp<-shet[is.element(shet$ensg, df2$ensg),]
shet_tmp_all<-shet[is.element(shet$ensg, df$ensg),]

A<-table(shet_tmp$shet_BIN)

random<-c()
for(I in 1:10000){
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
Program_P=min(P1,P2)

hoge<-data.frame(Program=Program, P=Program_P, meanG=MEANgamma_top100 - random_mean )
Program_P_sum=rbind(Program_P_sum, hoge)
}
Program_P_sum<-Program_P_sum[order(abs(Program_P_sum$meanG), decreasing=T),]
Program_P_sum<-Program_P_sum[order(Program_P_sum$P),]

Program_selected=paste0("P", Program_P_sum$Program[1:Program_N])

##step2. test regulator-burden correlation with step-wise regression to find the defined number of top enriched programs

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


##step3. evaluate the model prediction accuracy

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

df2<-merge(df[,"GENE", drop=F], GEP_reg_P, by="GENE")
row.names(df2)<-df2$GENE
df2<-df2[df$GENE,]

sum_effect<-data.frame(GENE=df$GENE, effect=0, regeffect=0, proeffect=0)
for(p in regulator_selected){
tmp<-res1[p, "Estimate"] *df[,p]
tmp[p.adjust(df2[,p], method="BH")>0.05]<-0
sum_effect$effect=sum_effect$effect + tmp
}

sum_effect$regeffect<-sign(sum_effect$effect)

for(p in rev(Program_selected)){  
tmp<-GEP[order(GEP[,p], decreasing=T), "GENE"][1:Program_top_def]

tmp2<-sign( Program_P_sum$meanG [Program_P_sum$Program==gsub("P", "", p)])
sum_effect$proeffect[is.element(sum_effect$GENE, tmp)]<-tmp2
}

LOF2<-data.frame(LOF, BIN=ifelse(abs(LOF$post_mean)>LOF_thresh, 1, 2))
for(BIN in 1) {
topgenes<-LOF2$GENE[LOF2$BIN == BIN]

tmp_mat<-LOF2[is.element(LOF2$GENE, topgenes), c("GENE", "post_mean")]

dff<-merge(tmp_mat, sum_effect, by="GENE")
dff<-data.frame(dff, post_mean_sign=sign(dff$post_mean))

dff<-data.frame(dff, predicted_sign=dff$proeffect)
dff$predicted_sign[dff$predicted_sign==0]<-dff$regeffect[dff$predicted_sign==0]


CONCO_N=length(which(dff$predicted_sign == dff$post_mean_sign))
DISCO_N=length(which(dff$predicted_sign == -dff$post_mean_sign & abs(dff$predicted_sign)>0 ))
TOTAL_N=nrow(dff)


##BG
tmp_mat2<-LOF2[!is.element(LOF2$GENE, topgenes), c("GENE", "post_mean")]

dff2<-merge(tmp_mat2, sum_effect, by="GENE")
dff2<-data.frame(dff2, post_mean_sign=sign(dff2$post_mean))

dff2<-data.frame(dff2, predicted_sign=dff2$proeffect)
dff2$predicted_sign[dff2$predicted_sign==0]<-dff2$regeffect[dff2$predicted_sign==0]

BG_CONCO_N=length(which(dff2$predicted_sign == dff2$post_mean_sign))
BG_DISCO_N=length(which(dff2$predicted_sign == -dff2$post_mean_sign & abs(dff2$predicted_sign)>0 ))
BG_TOTAL_N=nrow(dff2)

hoge<-data.frame(SEED="TRUE", BIN=BIN, CONCO_N=CONCO_N, DISCO_N=DISCO_N, TOTAL_N=TOTAL_N, BG_CONCO_N=BG_CONCO_N, BG_DISCO_N=BG_DISCO_N, BG_TOTAL_N=BG_TOTAL_N)
summary<-rbind(summary, hoge)
}

write.table(summary, paste0("permutation_test/P", Program_top_def, "/", TRAIT,  "_program", Program_N, "_regulator", Regulator_N, "_LOF", LOF_thresh, ".txt"), row.names=F, sep="\t", quote=F, append=F)


hoge<-data.frame(Program_selected=paste(Program_selected, collapse=","), Regulator_selected=paste(regulator_selected, collapse=","))
write.table(hoge, paste0("permutation_test/P", Program_top_def, "/SelectedPrograms_", TRAIT,  "_program", Program_N, "_regulator", Regulator_N, "_LOF", LOF_thresh, ".txt"), row.names=F, sep="\t", quote=F, append=F)


####
##permutation test
####

for(SEED in 1:20000){

##step1. test program burden effects to find the defined number of top enriched programs

if(TRAIT=="MCH"){
LOF=MCH
}
if(TRAIT=="RDW"){
LOF=RDW
}
if(TRAIT=="IRF"){
LOF=Ret
}

set.seed(SEED)
LOF$post_mean=sample(LOF$post_mean, nrow(LOF), replace=F)

Program_P_sum<-data.frame()
for(Program in 1:K){

tmp<-GEP[,c("GENE",paste0("P", Program))]
tmp<-tmp[order(tmp[,2], decreasing=T),]
colnames(tmp)[2]<-"GEPscore"
tmp<-data.frame(tmp, ensg=row.names(tmp))
df<-merge(tmp, LOF, by="ensg")

MEANgamma_top100=mean(df[order(abs(df$GEPscore), decreasing=T),"post_mean"][1:Program_top_def] )

##simplify to mwu test
Program_P=wilcox.test(df[order(abs(df$GEPscore), decreasing=T),"post_mean"][1:Program_top_def], df[order(abs(df$GEPscore), decreasing=T),"post_mean"][(Program_top_def + 1 ):nrow(df)])$p.value

hoge<-data.frame(Program=Program, P=Program_P, meanG=MEANgamma_top100 - mean(df$post_mean) )
Program_P_sum=rbind(Program_P_sum, hoge)
}
Program_P_sum<-Program_P_sum[order(abs(Program_P_sum$meanG), decreasing=T),]
Program_P_sum<-Program_P_sum[order(Program_P_sum$P),]

Program_selected=paste0("P", Program_P_sum$Program[1:Program_N])


##step2. regulator enrichment test to find the defined number of top enriched programs

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

##step3. calculate model enrichment

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

df2<-merge(df[,"GENE", drop=F], GEP_reg_P, by="GENE")
row.names(df2)<-df2$GENE
df2<-df2[df$GENE,]

sum_effect<-data.frame(GENE=df$GENE, effect=0, regeffect=0, proeffect=0)
for(p in regulator_selected){
tmp<-res1[p, "Estimate"] *df[,p]
tmp[p.adjust(df2[,p], method="BH")>0.05]<-0
sum_effect$effect=sum_effect$effect + tmp
}

sum_effect$regeffect<-sign(sum_effect$effect)

for(p in rev(Program_selected)){  
tmp<-GEP[order(GEP[,p], decreasing=T), "GENE"][1:Program_top_def]

tmp2<-sign( Program_P_sum$meanG [Program_P_sum$Program==gsub("P", "", p)])
sum_effect$proeffect[is.element(sum_effect$GENE, tmp)]<-tmp2
}

LOF2<-data.frame(LOF, BIN=ifelse(abs(LOF$post_mean)>LOF_thresh, 1, 2))
for(BIN in 1) {
topgenes<-LOF2$GENE[LOF2$BIN == BIN]

tmp_mat<-LOF2[is.element(LOF2$GENE, topgenes), c("GENE", "post_mean")]

dff<-merge(tmp_mat, sum_effect, by="GENE")
dff<-data.frame(dff, post_mean_sign=sign(dff$post_mean))

dff<-data.frame(dff, predicted_sign=dff$proeffect)
dff$predicted_sign[dff$predicted_sign==0]<-dff$regeffect[dff$predicted_sign==0]


CONCO_N=length(which(dff$predicted_sign == dff$post_mean_sign))
DISCO_N=length(which(dff$predicted_sign == -dff$post_mean_sign & abs(dff$predicted_sign)>0 ))
TOTAL_N=nrow(dff)


##BG
tmp_mat2<-LOF2[!is.element(LOF2$GENE, topgenes), c("GENE", "post_mean")]

dff2<-merge(tmp_mat2, sum_effect, by="GENE")
dff2<-data.frame(dff2, post_mean_sign=sign(dff2$post_mean))

dff2<-data.frame(dff2, predicted_sign=dff2$proeffect)
dff2$predicted_sign[dff2$predicted_sign==0]<-dff2$regeffect[dff2$predicted_sign==0]

BG_CONCO_N=length(which(dff2$predicted_sign == dff2$post_mean_sign))
BG_DISCO_N=length(which(dff2$predicted_sign == -dff2$post_mean_sign & abs(dff2$predicted_sign)>0 ))
BG_TOTAL_N=nrow(dff2)

hoge<-data.frame(SEED=SEED, BIN=BIN, CONCO_N=CONCO_N, DISCO_N=DISCO_N, TOTAL_N=TOTAL_N, BG_CONCO_N=BG_CONCO_N, BG_DISCO_N=BG_DISCO_N, BG_TOTAL_N=BG_TOTAL_N)
summary<-rbind(summary, hoge)

}
write.table(summary, paste0("permutation_test/P", Program_top_def, "/", TRAIT,  "_program", Program_N, "_regulator", Regulator_N, "_LOF", LOF_thresh, ".txt"), row.names=F, sep="\t", quote=F, append=F)

}

write.table(summary, paste0("permutation_test/P", Program_top_def, "/", TRAIT,  "_program", Program_N, "_regulator", Regulator_N, "_LOF", LOF_thresh, ".txt"), row.names=F, sep="\t", quote=F, append=F)

