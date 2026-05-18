args <- commandArgs(trailingOnly = TRUE)
burden_file <- if (length(args) >= 1) as.character(args[1]) else "Backman_2021_86_M1_001.summary_statistics.csv"

LOF<-read.table(paste0("data/LoF/raw_burden/", burden_file), header=T, stringsAsFactor=F, sep=",")
LOF<-data.frame(LOF, Z=LOF$beta/LOF$standard_error, P=2*pnorm(-abs(LOF$beta/LOF$standard_error)))

P<-seq(1/(nrow(LOF)+1), 1-1/(nrow(LOF)+1), by=1/(nrow(LOF)+1))

df<-LOF[,c("beta", "P", "ensg")]
df<-data.frame(df, FDR=p.adjust(df$P, method="BH"))
corresp<-read.table("data/gencode_v41_gname_gid_ALL_sorted_onlyID", header=F, stringsAsFactor=F, row.names=1)
df<-data.frame(df, LABEL=corresp[df$ensg, 1], geneset="other")

for(geneset in rev(c("HALLMARK_HEME_METABOLISM", "Hematopoiesisgenes", "mitotic_cell_cycle", "positive_macromolecule_synthesis" ))){
heme<-read.table(paste0("data/geneset/", geneset, ".txt"), header=F, stringsAsFactor=F)[,1]
df$geneset[is.element(df$LABEL, heme)]<-geneset
}

df$LABEL[abs(df$FDR)>0.1]<-""
df$LABEL[df$geneset=="other"]<-""

df_heme<-df[df$geneset!="other" & abs(df$FDR)<0.1,]
df_other<-df[df$geneset=="other",]

library(ggplot2)
library(ggrepel)
library(ggsci)
library(ggbreak) 

df$geneset<-gsub("HALLMARK_HEME_METABOLISM", "Heme metabolism", df$geneset)
df$geneset<-gsub("Hematopoiesisgenes", "Hematopoiesis", df$geneset)
df$geneset<-gsub("mitotic_cell_cycle", "Mitotic cell cycle", df$geneset)
df$geneset<-gsub("positive_macromolecule_synthesis", "Macromolecule byosynthesis, positive reg", df$geneset)

df_heme$geneset<-gsub("HALLMARK_HEME_METABOLISM", "Heme metabolism", df_heme$geneset)
df_heme$geneset<-gsub("Hematopoiesisgenes", "Hematopoiesis", df_heme$geneset)
df_heme$geneset<-gsub("mitotic_cell_cycle", "Mitotic cell cycle",  df_heme$geneset)
df_heme$geneset<-gsub("positive_macromolecule_synthesis", "Macromolecule byosynthesis, positive reg", df_heme$geneset)

thresh<-sort(df$P[df$FDR<0.1], decreasing=T)[1]

g<-ggplot(df, aes(x=beta, y=-log10(P), label=LABEL,color=geneset))
g<-g +  theme_classic(base_size = 40, base_family = "Helvetica")
g<-g + geom_point(size=3, alpha=0.5)
g<-g + geom_point(data=df_heme, size=3, alpha=0.8)
g<-g + geom_hline(yintercept=-log10(thresh), linetype="dashed")
g<-g + geom_text_repel(data=subset(df_heme, beta<0 & FDR<0.01), max.overlaps=200,  point.size=3, size=10, show.legend = FALSE, xlim=c(-Inf, -0.6) )
g<-g + geom_text_repel(data=subset(df_heme, beta>0 & FDR<0.01), max.overlaps=200,  point.size=3, size=10, show.legend = FALSE, xlim=c(0.6, Inf) )
g<-g + scale_color_manual(values = c("other"="grey60", "Heme metabolism"=pal_d3("category10")(4)[1], "Hematopoiesis"=pal_d3("category10")(4)[2],
  "Mitotic cell cycle"=pal_d3("category10")(4)[3], "Macromolecule byosynthesis, positive reg"=pal_d3("category10")(4)[4]),
breaks=c("Heme metabolism", "Hematopoiesis", "Mitotic cell cycle", "Macromolecule byosynthesis, positive reg", "other" ))
g<-g + ylab("-log10 (P)") +  xlab("Effect size")
g<-g + scale_x_continuous(limits=c(-3,2.5), sec.axis = dup_axis(breaks = NULL))
g<-g + scale_y_break(c(60, 220), scales = 0.05)
ggsave(plot=g, "Fig2B.pdf", width=20, height=23,  dpi=600)
