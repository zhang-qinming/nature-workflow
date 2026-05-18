args   <- commandArgs(trailingOnly = T)  
options("scipen"=10)

K<-as.numeric( args[1] ) ##60

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


##MCH plots

Program=16


tmpA<-GEP_reg_beta$GENE[GEP_reg_beta[,paste0("P", Program)]>0]
tmpB<-GEP_reg_P$GENE[ p.adjust(GEP_reg_P[,paste0("P", Program)], method="BH") <0.1]

Perturb_pos_genes<-intersect(tmpA, tmpB)

tmpA<-GEP_reg_beta$GENE[GEP_reg_beta[,paste0("P", Program)]<=0]
tmpB<-GEP_reg_P$GENE[ p.adjust(GEP_reg_P[,paste0("P", Program)], method="BH") <0.1]

Perturb_neg_genes<-intersect(tmpA, tmpB)

nonsig_genes<-GEP_reg_P$GENE[ p.adjust(GEP_reg_P[,paste0("P", Program)], method="BH") > 0.1 ]

df_MCH1<-MCH[is.element(MCH$gene, Perturb_pos_genes),  c("gene", "post_mean")]
df_MCH2<-MCH[is.element(MCH$gene, Perturb_neg_genes),  c("gene", "post_mean")]
df_MCH3<-MCH[is.element(MCH$gene, nonsig_genes),  c("gene", "post_mean")]


tmp.random1<-data.frame()
for(i in 1:5000){
random<-sort(sample(df_MCH3$post_mean, nrow(df_MCH1), replace=F))
tmp.random1<-rbind(tmp.random1, t(data.frame(num=sort(random))))
}
expected=apply(tmp.random1, 2, median)
low_90=apply(tmp.random1, 2, function(x){sort(x)[250]})
high_90=apply(tmp.random1, 2, function(x){sort(x)[4750]})
df.random1=data.frame(expected=expected, low_90=low_90, high_90=high_90, LABEL="")

tmp.random2<-data.frame()
for(i in 1:5000){
random<-sort(sample(df_MCH3$post_mean, nrow(df_MCH2), replace=F))
tmp.random2<-rbind(tmp.random2, t(data.frame(num=sort(random))))
}
expected=apply(tmp.random2, 2, median)
low_90=apply(tmp.random2, 2, function(x){sort(x)[250]})
high_90=apply(tmp.random2, 2, function(x){sort(x)[4750]})
df.random2=data.frame(expected=expected, low_90=low_90, high_90=high_90, LABEL="")


df_MCH1<-df_MCH1[order(df_MCH1$post_mean),]
df_MCH1<-data.frame(df_MCH1, expected=df.random1$expected, low_90=df.random1$low_90, high_90=df.random1$high_90)

df_MCH2<-df_MCH2[order(df_MCH2$post_mean),]
df_MCH2<-data.frame(df_MCH2, expected=df.random2$expected, low_90=df.random2$low_90, high_90=df.random2$high_90)


df_MCH1<-data.frame(df_MCH1, LABEL=df_MCH1$gene)
df_MCH1$LABEL[abs(df_MCH1$post_mean)<0.15]<-""
df_MCH1$LABEL[df_MCH1$post_mean>0]<-""
df_MCH12<-df_MCH1[abs(df_MCH1$post_mean)>0.15 & df_MCH1$post_mean<0,]

df_MCH2<-data.frame(df_MCH2, LABEL=df_MCH2$gene)
df_MCH2$LABEL[abs(df_MCH2$post_mean)<0.15]<-""
df_MCH2$LABEL[df_MCH2$post_mean<0]<-""
df_MCH22<-df_MCH2[!is.element(df_MCH2$LABEL, ""),]


library(ggplot2)
library(ggrepel)
library(cowplot)

g1<-ggplot(df_MCH2, aes(x=expected, y=post_mean, label=LABEL))
g1<-g1 +  theme_classic(base_size = 40, base_family = "Helvetica")
g1<-g1+ geom_point(color="black", fill="grey60", shape=21, alpha=0.8, size=4)
g1<-g1+ geom_point(data=df_MCH22, color="black", fill="grey60", shape=21, alpha=0.8, size=4)


g2<-ggplot(df_MCH1, aes(x=expected, y=post_mean, label=LABEL))
g2<-g2 +  theme_classic(base_size = 40, base_family = "Helvetica")
g2<-g2+ geom_point(color="black", fill="grey60", shape=21, alpha=0.8, size=4)
g2<-g2+ geom_point(data=df_MCH12,color="black", fill="grey60", shape=21, alpha=0.8, size=4)

g1<-g1 + geom_text_repel(nudge_x=0.25, size=8)
g2<-g2 + geom_text_repel( max.overlaps=100, size=8, box.padding=0.5)

g1<-g1 + geom_abline(intercept=0, slope=1, linetype="dashed")
g2<-g2 + geom_abline(intercept=0, slope=1, linetype="dashed")

g1<-g1 + ylab(expression(paste("O[Burden test ", gamma, ", MCH]"))) + xlab(expression(paste("E[Burden test ", gamma, ", MCH]")))
g2<-g2 + ylab(expression(paste("O[Burden test ", gamma, ", MCH]"))) + xlab(expression(paste("E[Burden test ", gamma, ", MCH]")))
g3<-plot_grid(g2, g1, nrow=1)
ggsave(plot=g2, paste0("DownRegulatorgenes_", Program, "_MCH_Fig5J.pdf"), width=8, height=10)
ggsave(plot=g1, paste0("UpRegulatorgenes_", Program, "_MCH_Fig5J.pdf"), width=8, height=10)



####
##highlight autophagosome genes in the program genes plot
##GO:0005776 autophagosome

geneset<-read.table("data/geneset/Autophagosome_genes.txt", header=F, stringsAsFactor=F)[,1]

NUM=200

tmp<-GEP[,c("GENE",paste0("P", Program))]
tmp<-tmp[order(tmp[,2], decreasing=T),]
genes<-tmp$GENE[1:NUM]
genes<-na.omit(genes)

BG<-intersect(MCH$gene, corresp[is.element(corresp[,2], GEP$GENE),2])
df<-MCH[is.element(MCH$gene, BG),]
tmp<-data.frame(tmp=ifelse(is.element(df$gene, genes), 1,0))
colnames(tmp)<-paste0("GEP")

df<-cbind(df, tmp)
df[,paste0("GEP")]<-ordered(df[,paste0("GEP")])
library(MASS)

BG<-BG[!is.element(BG, genes)]
df_MCH1<-MCH[is.element(MCH$gene, genes),  c("gene", "post_mean")]
df_MCH3<-MCH[is.element(MCH$gene, BG),  c("gene", "post_mean")]

df_MCH1<-df_MCH1[order(df_MCH1$post_mean),]

tmp.random<-data.frame()
for(i in 1:5000){
random<-sort(sample(df_MCH3$post_mean, nrow(df_MCH1), replace=F))
tmp.random<-rbind(tmp.random, t(data.frame(num=sort(random))))
}
expected=apply(tmp.random, 2, median)
low_95=apply(tmp.random, 2, function(x){sort(x)[250]})
high_95=apply(tmp.random, 2, function(x){sort(x)[9750]})

df.random=data.frame(expected=expected, low_95=low_95, high_95=high_95, LABEL="")
expected=sort(apply(tmp.random, 2, median))
df_MCH1<-data.frame(df_MCH1, expected=expected)

df_MCH1<-data.frame(df_MCH1, LABEL=df_MCH1$gene)
df_MCH1$LABEL[!is.element(df_MCH1$gene, geneset)]<-""
df_MCH12<-df_MCH1[!is.element(df_MCH1$LABEL, ""),]


library(ggplot2)
library(ggrepel)
g<-ggplot(df_MCH1, aes(x=expected, y=post_mean, label=LABEL))
g<-g +  theme_classic(base_size = 40, base_family = "Helvetica")
g<-g+ geom_point(color="black", fill="grey60", shape=21, alpha=0.8, size=4)
g<-g+ geom_point(data=df_MCH12, color="black", fill="grey60", shape=21, alpha=0.8, size=4)
g<-g + ylab("gene effect size")
g<-g + geom_text_repel(size=7, max.overlaps=1000,  box.padding=0.5)
g<-g + geom_abline(intercept=0, slope=1, linetype="dashed")
g<-g + ylab(expression(paste("O[Burden test ", gamma, ", MCH]"))) + xlab(expression(paste("E[Burden test ", gamma, ", MCH]")))
ggsave(plot=g, paste0("Programgenes_", Program, "_MCH_Fig5J.pdf"), width=8, height=10)



##RDW plots

Program=16


tmpA<-GEP_reg_beta$GENE[GEP_reg_beta[,paste0("P", Program)]>0]
tmpB<-GEP_reg_P$GENE[ p.adjust(GEP_reg_P[,paste0("P", Program)], method="BH") <0.1]

Perturb_pos_genes<-intersect(tmpA, tmpB)

tmpA<-GEP_reg_beta$GENE[GEP_reg_beta[,paste0("P", Program)]<=0]
tmpB<-GEP_reg_P$GENE[ p.adjust(GEP_reg_P[,paste0("P", Program)], method="BH") <0.1]

Perturb_neg_genes<-intersect(tmpA, tmpB)

nonsig_genes<-GEP_reg_P$GENE[ p.adjust(GEP_reg_P[,paste0("P", Program)], method="BH") > 0.1 ]

df_RDW1<-RDW[is.element(RDW$gene, Perturb_pos_genes),  c("gene", "post_mean")]
df_RDW2<-RDW[is.element(RDW$gene, Perturb_neg_genes),  c("gene", "post_mean")]
df_RDW3<-RDW[is.element(RDW$gene, nonsig_genes),  c("gene", "post_mean")]


tmp.random1<-data.frame()
for(i in 1:5000){
random<-sort(sample(df_RDW3$post_mean, nrow(df_RDW1), replace=F))
tmp.random1<-rbind(tmp.random1, t(data.frame(num=sort(random))))
}
expected=apply(tmp.random1, 2, median)
low_90=apply(tmp.random1, 2, function(x){sort(x)[250]})
high_90=apply(tmp.random1, 2, function(x){sort(x)[4750]})
df.random1=data.frame(expected=expected, low_90=low_90, high_90=high_90, LABEL="")

tmp.random2<-data.frame()
for(i in 1:5000){
random<-sort(sample(df_RDW3$post_mean, nrow(df_RDW2), replace=F))
tmp.random2<-rbind(tmp.random2, t(data.frame(num=sort(random))))
}
expected=apply(tmp.random2, 2, median)
low_90=apply(tmp.random2, 2, function(x){sort(x)[250]})
high_90=apply(tmp.random2, 2, function(x){sort(x)[4750]})
df.random2=data.frame(expected=expected, low_90=low_90, high_90=high_90, LABEL="")


df_RDW1<-df_RDW1[order(df_RDW1$post_mean),]
df_RDW1<-data.frame(df_RDW1, expected=df.random1$expected, low_90=df.random1$low_90, high_90=df.random1$high_90)

df_RDW2<-df_RDW2[order(df_RDW2$post_mean),]
df_RDW2<-data.frame(df_RDW2, expected=df.random2$expected, low_90=df.random2$low_90, high_90=df.random2$high_90)


df_RDW1<-data.frame(df_RDW1, LABEL=df_RDW1$gene)
df_RDW1$LABEL[abs(df_RDW1$post_mean)<0.2]<-""
df_RDW1$LABEL[df_RDW1$post_mean<0]<-""
df_RDW12<-df_RDW1[abs(df_RDW1$post_mean)>0.25 & df_RDW1$post_mean>0,]

df_RDW2<-data.frame(df_RDW2, LABEL=df_RDW2$gene)
df_RDW2$LABEL[abs(df_RDW2$post_mean)<0.12]<-""
df_RDW2$LABEL[df_RDW2$post_mean>0]<-""
df_RDW22<-df_RDW2[!is.element(df_RDW2$LABEL, ""),]


library(ggplot2)
library(ggrepel)
library(cowplot)

g1<-ggplot(df_RDW2, aes(x=expected, y=post_mean, label=LABEL))
g1<-g1 +  theme_classic(base_size = 40, base_family = "Helvetica")
g1<-g1+ geom_point(color="black", fill="grey60", shape=21, alpha=0.8, size=4)
g1<-g1+ geom_point(data=df_RDW22, color="black", fill="grey60", shape=21, alpha=0.8, size=4)


g2<-ggplot(df_RDW1, aes(x=expected, y=post_mean, label=LABEL))
g2<-g2 +  theme_classic(base_size = 40, base_family = "Helvetica")
g2<-g2+ geom_point(color="black", fill="grey60", shape=21, alpha=0.8, size=4)
g2<-g2+ geom_point(data=df_RDW12,color="black", fill="grey60", shape=21, alpha=0.8, size=4)

g1<-g1 + geom_text_repel(nudge_x=0.25, size=8)
g2<-g2 + geom_text_repel( max.overlaps=100, size=8, box.padding=0.5)

g1<-g1 + geom_abline(intercept=0, slope=1, linetype="dashed")
g2<-g2 + geom_abline(intercept=0, slope=1, linetype="dashed")

g1<-g1 + ylab(expression(paste("O[Burden test ", gamma, ", RDW]"))) + xlab(expression(paste("E[Burden test ", gamma, ", RDW]")))
g2<-g2 + ylab(expression(paste("O[Burden test ", gamma, ", RDW]"))) + xlab(expression(paste("E[Burden test ", gamma, ", RDW]")))
g3<-plot_grid(g2, g1, nrow=1)
ggsave(plot=g2, paste0("DownRegulatorgenes_", Program, "_RDW_Fig5J.pdf"), width=8, height=10)
ggsave(plot=g1, paste0("UpRegulatorgenes_", Program, "_RDW_Fig5J.pdf"), width=8, height=10)



####
##highlight autophagosome genes in the program genes plot
##GO:0005776 autophagosome

geneset<-read.table("data/geneset/Autophagosome_genes.txt", header=F, stringsAsFactor=F)[,1]

NUM=200

tmp<-GEP[,c("GENE",paste0("P", Program))]
tmp<-tmp[order(tmp[,2], decreasing=T),]
genes<-tmp$GENE[1:NUM]
genes<-na.omit(genes)

BG<-intersect(RDW$gene, corresp[is.element(corresp[,2], GEP$GENE),2])
df<-RDW[is.element(RDW$gene, BG),]
tmp<-data.frame(tmp=ifelse(is.element(df$gene, genes), 1,0))
colnames(tmp)<-paste0("GEP")

df<-cbind(df, tmp)
df[,paste0("GEP")]<-ordered(df[,paste0("GEP")])
library(MASS)

BG<-BG[!is.element(BG, genes)]
df_RDW1<-RDW[is.element(RDW$gene, genes),  c("gene", "post_mean")]
df_RDW3<-RDW[is.element(RDW$gene, BG),  c("gene", "post_mean")]

df_RDW1<-df_RDW1[order(df_RDW1$post_mean),]

tmp.random<-data.frame()
for(i in 1:5000){
random<-sort(sample(df_RDW3$post_mean, nrow(df_RDW1), replace=F))
tmp.random<-rbind(tmp.random, t(data.frame(num=sort(random))))
}
expected=apply(tmp.random, 2, median)
low_95=apply(tmp.random, 2, function(x){sort(x)[250]})
high_95=apply(tmp.random, 2, function(x){sort(x)[9750]})

df.random=data.frame(expected=expected, low_95=low_95, high_95=high_95, LABEL="")
expected=sort(apply(tmp.random, 2, median))
df_RDW1<-data.frame(df_RDW1, expected=expected)

df_RDW1<-data.frame(df_RDW1, LABEL=df_RDW1$gene)
df_RDW1$LABEL[!is.element(df_RDW1$gene, geneset)]<-""
df_RDW12<-df_RDW1[!is.element(df_RDW1$LABEL, ""),]


library(ggplot2)
library(ggrepel)
g<-ggplot(df_RDW1, aes(x=expected, y=post_mean, label=LABEL))
g<-g +  theme_classic(base_size = 40, base_family = "Helvetica")
g<-g+ geom_point(color="black", fill="grey60", shape=21, alpha=0.8, size=4)
g<-g+ geom_point(data=df_RDW12, color="black", fill="grey60", shape=21, alpha=0.8, size=4)
g<-g + ylab("gene effect size")
g<-g + geom_text_repel(size=7, max.overlaps=1000,  box.padding=0.5)
g<-g + geom_abline(intercept=0, slope=1, linetype="dashed")
g<-g + ylab(expression(paste("O[Burden test ", gamma, ", RDW]"))) + xlab(expression(paste("E[Burden test ", gamma, ", RDW]")))
ggsave(plot=g, paste0("Programgenes_", Program, "_RDW_Fig5J.pdf"), width=8, height=10)

