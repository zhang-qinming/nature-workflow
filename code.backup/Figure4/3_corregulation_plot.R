args   <- commandArgs(trailingOnly = T)  
options("scipen"=10)

T1<-as.character( args[1] ) ##e.g., "P25"
T2<-as.character( args[2] ) ##e.g., "P16"

K<-as.numeric( args[3] ) ##60

##data loading

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
corresp<-read.table("data/gencode_v41_gname_gid_ALL_sorted_onlyID", header=F, stringsAsFactor=F)
corresp<-corresp[!duplicated(corresp[,1]),]
row.names(corresp)<-corresp[,1]

df<-GEP_reg_beta

library(ggplot2)
library(ggallin)
target<-unique(c(GEP_reg_P$GENE[p.adjust(GEP_reg_P[,T1], method="BH")<0.05]))

df<-df[is.element(df$GENE, target),]

df_f<-df[,c(T1,T2, "GENE")]
colnames(df_f)[1:2]<-c("x", "y")

## df_f$y[df_f$y>3]<-3  ##if needed, truncate the axis.

g<-ggplot(df_f, aes(x=x, y=y))
g<-g +  theme_classic(base_size = 30, base_family = "Helvetica")
g<-g +  geom_point(shape=21, color="black", fill=NA,stroke=1, size=3, alpha=0.8)
g<-g + geom_vline(xintercept=0, linetype="dashed")
g<-g + geom_hline(yintercept=0, linetype="dashed")
g<-g + xlab(paste0("effect size on ", T1)) + ylab(paste0("effect size on ", T2))
g<-g + geom_smooth(method = "loess",method.args= list(degree = 1), span = 0.25)
g<-g + theme(legend.position="none")
ggsave(plot=g, paste0( T1, "_vs_", T2, "_Fig5.pdf"), width=8, height=8 )

