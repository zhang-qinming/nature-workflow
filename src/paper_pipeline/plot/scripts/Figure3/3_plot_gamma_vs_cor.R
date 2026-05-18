args <- commandArgs(trailingOnly = TRUE)
FILE <- if (length(args) >= 1) as.character(args[1]) else "Backman_2021_86.per_gene_estimates.tsv"

data<-read.table(paste0("data/Perturbseq/trait_association/K562GW/GeneLevel/", FILE, "_geneRegulation_correlation.txt"), header=T, stringsAsFactor=F)
LOF<-read.table(paste0("data/LoF/GeneBayes_posterior/", FILE), header=T, stringsAsFactor=F)
df<-merge(data, LOF, by="ensg")
df<-data.frame(df, LABEL=df$gene)
df<-df[order(abs(df$post_mean), decreasing=T),]
df$LABEL[-grep("^HB", df$LABEL)]<-""
df$LABEL[!is.element(df$LABEL, c("HBA1"))]<-""

df<-data.frame(df, signedP=sign(df$beta_withShet)*(-log10(df$P_withShet)))
df<-df[order(df$P_withShet, decreasing=T),]
THRESH=-log10(df$P_withShet[ p.adjust(df$P_withShet, method="BH")<0.05 ] [1])
df_HBA1<-df[is.element(df$gene, "HBA1"),,drop=F]
df$LABEL<-gsub("HBA1", "", df$LABEL)

library(ggplot2)
library(ggrepel)
g<-ggplot(df, aes(x=post_mean, y=signedP, label=LABEL))
g<-g +  theme_classic(base_size = 30, base_family = "Helvetica")
g<-g + geom_point(size=3)
g<-g + geom_point(data=df_HBA1, size=5, color="red")
g<-g + geom_text_repel(data=df_HBA1, size=10, point.size=5)
g<-g + ylab("Burden-regulator correlation, signed -log10(P)") + xlab(paste0("Burden test gamma for ", FILE))
g<-g + theme(legend.position = "none")
g<-g + geom_hline(yintercept=THRESH, linetype="dashed")
g<-g + geom_hline(yintercept=-THRESH, linetype="dashed")
g<-g + scale_y_continuous(breaks=c(-8,-4,0,4,8), limits=c(-8,8))
ggsave(plot=g, "Fig3c.pdf", width=12, height=12)
