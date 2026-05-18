args   <- commandArgs(trailingOnly = T) 

trait_name <- as.character( args[1] ) ##"MCH"

##GWAS all top variants, closest genes
GWAS_go<-read.table(paste0("enrichment_result/", trait_name, "_",  "GWAS_GO.txt"), header=T, sep="\t", quote="", stringsAsFactor=F)
GWAS_msig<-read.table(paste0("enrichment_result/", trait_name, "_",  "GWAS_MsigDb.txt"), header=T, sep="\t", quote="",stringsAsFactor=F)

##LOF, after GeneBayes
LOF_go<-read.table( paste0("enrichment_result/", trait_name, "_",  "LOF_GB_GO.txt"), header=T,sep="\t", quote="", stringsAsFactor=F)
LOF_msig<-read.table( paste0("enrichment_result/", trait_name, "_",  "LOF_GB_MsigDb.txt"), header=T,sep="\t", quote="", stringsAsFactor=F)

##

GWAS_df<-rbind(GWAS_go[,c("ID", "Description", "pvalue")], GWAS_msig[,c("ID", "Description", "pvalue")])
LOF_df<-rbind(LOF_go[,c("ID", "Description", "pvalue")], LOF_msig[,c("ID", "Description", "pvalue")])

colnames(GWAS_df)[3]<-"GWAS_pvalue"
colnames(LOF_df)[3]<-"LOF_pvalue"

df<-merge(GWAS_df, LOF_df, by="ID")

df<-data.frame(df, LABEL=df$Description.x)
df$LABEL<-gsub("HALLMARK_HEME_METABOLISM", "heme metabolism", df$LABEL)

##change colors
library(ggplot2)
library(ggrepel)
library(ggsci)

df$LABEL[!is.element(df$LABEL, c("heme metabolism", "mitotic cell cycle", "hematopoietic or lymphoid organ development", "hemopoiesis",
  "positive regulation of macromolecule biosynthetic process", "nuclear protein-containing complex", "autophagy" ))]<-""
df$LABEL<-gsub("positive regulation of macromolecule biosynthetic process", "Macromolecule byosynthesis, positive reg", df$LABEL)
df$LABEL<-gsub("heme metabolism", "Heme metabolism", df$LABEL)
df$LABEL<-gsub("hemopoiesis", "Hematopoiesis", df$LABEL)
df$LABEL<-gsub("mitotic cell cycle", "Mitotic cell cycle", df$LABEL)

df$COLOR=df$LABEL
df$COLOR[is.element(df$COLOR, "")]<-"other"
df$COLOR[!is.element(df$COLOR, c("Heme metabolism", "Mitotic cell cycle", "Hematopoiesis",
  "Macromolecule byosynthesis, positive reg", "other" ))]<-"other2"

df2<-df[df$LABEL != "",]

g<-ggplot(df, aes(x=-log10(GWAS_pvalue), y=-log10(LOF_pvalue), label=LABEL, color=COLOR))
g<-g +  theme_classic(base_size = 50, base_family = "Helvetica")
g<-g + geom_point(alpha=1, size=6)
g<-g + xlab("GWAS enrichment, -log10(P)") + ylab("LOF enrichment, -log10(P)")
g<-g + geom_point(data=df2, alpha=1, size=6)
g<-g + geom_text_repel(data=df, size=8, point.size=6)
g<-g + scale_color_manual(values = c("other"="grey60", "other2"="black", "Heme metabolism"=pal_d3("category10")(4)[1], "Hematopoiesis"=pal_d3("category10")(4)[2],
  "Mitotic cell cycle"=pal_d3("category10")(4)[3], "Macromolecule byosynthesis, positive reg"=pal_d3("category10")(4)[4]),  guide="none")
ggsave(plot=g, "Fig2C.pdf", width=14, height=14)
