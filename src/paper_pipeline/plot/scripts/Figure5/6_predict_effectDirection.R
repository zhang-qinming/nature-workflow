args   <- commandArgs(trailingOnly = T)   # enable reading factors from linux command line


TRAIT=as.character( args[1] ) #"MCH", "RDW", "IRF"
Program_top_def=as.numeric( args[2] ) ##100,200

K=if (length(args) >= 3) as.numeric(args[3]) else 60
mch_posterior_path <- if (length(args) >= 4) as.character(args[4]) else "data/LoF/GeneBayes_posterior/Backman_2021_86.per_gene_estimates.tsv"
rdw_posterior_path <- if (length(args) >= 5) as.character(args[5]) else "data/LoF/GeneBayes_posterior/Backman_2021_88.per_gene_estimates.tsv"
irf_posterior_path <- if (length(args) >= 6) as.character(args[6]) else "data/LoF/GeneBayes_posterior/Backman_2021_106.per_gene_estimates.tsv"
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

##step1. calculate program mean gamma


Program_P_sum<-data.frame()
for(Program in 1:K){

tmp<-GEP[,c("GENE",paste0("P", Program))]
tmp<-tmp[order(tmp[,2], decreasing=T),]
colnames(tmp)[2]<-"GEPscore"
tmp<-data.frame(tmp, ensg=row.names(tmp))
df<-merge(tmp, LOF, by="ensg")
df2<-df[order(df$GEPscore, decreasing=T),][1:Program_top_def,]

MEANgamma_top100=mean(df2[order(abs(df2$post_mean), decreasing=T),"post_mean"][1:Program_top_def] )

 
##make random 10000 set of genes, matched to shet bin

shet_tmp<-shet[is.element(shet$ensg, df2$ensg),]
shet_tmp_all<-shet[is.element(shet$ensg, df$ensg),]

A<-table(shet_tmp$shet_BIN)

random<-c()
for(I in 1:1000){
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

hoge<-data.frame(Program=Program, P=program_P, meanG=MEANgamma_top100 - random_mean )
Program_P_sum=rbind(Program_P_sum, hoge)
}
Program_P_sum<-Program_P_sum[order(abs(Program_P_sum$meanG), decreasing=T),]
Program_P_sum<-Program_P_sum[order(Program_P_sum$P),]



##step3. calculate model enrichment
if(TRAIT=="MCH"){
regulator_selected=c("P25", "P40", "P16")
Program_selected=c("P40", "P25", "P4", "P13", "P17")
}

if(TRAIT=="RDW"){
regulator_selected=c("P25", "P43", "P16")
Program_selected=c("P40", "P1")
}

if(TRAIT=="IRF"){
regulator_selected=c("P43")
Program_selected=c("P43", "P5")
}


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

tmp3<-tmp[!is.element(tmp, sum_effect$GENE)]
if(length(tmp3)>0){
sum_effect<-rbind(sum_effect, data.frame(GENE=tmp3, effect=0, regeffect=0, proeffect=tmp2))
}}


LOF2<-data.frame(LOF, BIN=ifelse(abs(LOF$post_mean)>0.1, 1, ifelse(abs(LOF$post_mean)>0.08, 2, ifelse(abs(LOF$post_mean)>0.054, 3, 
	ifelse(abs(LOF$post_mean)>0.05, 4, 5)))))

BIN=1

topgenes<-LOF2$GENE[LOF2$BIN == BIN]

tmp_mat<-LOF2[is.element(LOF2$GENE, topgenes), c("GENE", "post_mean")]

dff<-merge(tmp_mat, sum_effect, by="GENE")
dff<-data.frame(dff, post_mean_sign=sign(dff$post_mean))

dff<-data.frame(dff, predicted_sign=dff$proeffect)
dff$predicted_sign[dff$predicted_sign==0]<-dff$regeffect[dff$predicted_sign==0]


CONCO_N=length(which(dff$predicted_sign == dff$post_mean_sign))
DISCO_N=length(which(dff$predicted_sign == -dff$post_mean_sign & abs(dff$predicted_sign)>0 ))
TOTAL_N=nrow(dff)

mat<-dff

##annotate with regulatory effects and program genes

for(p in c(regulator_selected)){
target<-intersect(GEP_reg_P$GENE[p.adjust(GEP_reg_P[,p], method="BH")<0.05], mat$GENE)
mat<-data.frame(mat, X=0)
colnames(mat)[ncol(mat)]<-paste0(p, "reg")
for(g in target){
mat[is.element(mat$GENE, g),ncol(mat)] <- sign(GEP_reg_beta[,p][is.element(GEP_reg_beta$GENE, g)])
}
}

for(p in Program_selected){
target<-intersect(GEP[order(GEP[,p], decreasing=T), "GENE"][1:Program_top_def], mat$GENE)
tmp2<-sign( Program_P_sum$meanG [Program_P_sum$Program==gsub("P", "", p)])
mat<-data.frame(mat, X=0)
colnames(mat)[ncol(mat)]<-paste0(p, "pro")
for(g in target){
mat[is.element(mat$GENE, g),ncol(mat)] <- tmp2
}
}

write.table(mat, paste0(TRAIT, "_all_gamma0.1_genes_prediction.txt"), row.names=F, sep="\t", quote=F)

mat2<-mat[mat$predicted_sign == mat$post_mean_sign,]
write.table(mat2, paste0(TRAIT, "_concordant_gamma0.1_genes_prediction.txt"), row.names=F, sep="\t", quote=F)

##use this for making gene-to-program-to-trait map
