args <- commandArgs(trailingOnly = TRUE)
K <- if (length(args) >= 1) as.numeric(args[1]) else 60

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

##
mat<-combn(K,2)
summary<-data.frame()
for(p in 1:ncol(mat)){
P_A=paste0("P", mat[1,p])
P_B=paste0("P", mat[2,p])
library(DescTools)
genesA<-GEP_reg_P$GENE[p.adjust(GEP_reg_P[,P_A], method="BH")<0.05]
Rhox= cor.test(GEP_reg_beta[is.element(GEP_reg_beta$GENE, genesA), P_A], GEP_reg_beta[is.element(GEP_reg_beta$GENE, genesA), P_B], method="spearman")$estimate

Zx=FisherZ(Rhox)
Nx=length(genesA)

genesB<-GEP_reg_P$GENE[p.adjust(GEP_reg_P[,P_B], method="BH")<0.05]
Rhoy= cor.test(GEP_reg_beta[is.element(GEP_reg_beta$GENE, genesB), P_A], GEP_reg_beta[is.element(GEP_reg_beta$GENE, genesB), P_B], method="spearman")$estimate

Zy=FisherZ(Rhoy)
Ny=length(genesB)

aic_1 <- 2*1 - (dnorm(Zx, mean=Zx, sd=sqrt(1/(Nx-3)), log=T) + dnorm(Zy, mean=0, sd=sqrt(1/(Ny-3)), log=T))
aic_2 <- 2*1 - (dnorm(Zx, mean=0, sd=sqrt(1/(Nx-3)), log=T) + dnorm(Zy, mean=Zy, sd=sqrt(1/(Ny-3)), log=T))
aic_3 <- 2*0 - (dnorm(Zx, mean=0, sd=sqrt(1/(Nx-3)), log=T) + dnorm(Zy, mean=0, sd=sqrt(1/(Ny-3)), log=T))
weighted_mean_z=(Zx * (Nx-3) + Zy *(Ny-3)) / (Nx-3 + Ny-3)
aic_4 <- 2*1 - (dnorm(Zx, mean=weighted_mean_z, sd=sqrt(1/(Nx-3)), log=T) + dnorm(Zy, mean=weighted_mean_z, sd=sqrt(1/(Ny-3)), log=T))

aic_causal<-min(aic_1, aic_2)
aic_noncausal<-min(aic_3, aic_4)

r=exp((aic_causal - aic_noncausal)/2)
##defined <0.05 as causal

hoge<-data.frame(P_A=P_A, P_B=P_B, Rhox=Rhox, Rhoy=Rhoy, aic_causal=aic_causal, aic_noncausal=aic_noncausal, likelihood=r)
summary<-rbind(summary, hoge)
}

write.table(summary, "data/programs_causal_relation_summary.txt", row.names=F, sep="\t", quote=F)
