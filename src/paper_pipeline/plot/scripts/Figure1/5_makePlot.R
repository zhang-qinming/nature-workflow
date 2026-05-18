args <- commandArgs(trailingOnly = TRUE)

##1. p values for enrichment to K562, RBC traits colored.
wd <- if (length(args) >= 1) as.character(args[1]) else "data/LDSC"
setwd(wd)
dir.create("figure")

data<-read.table("result_summary/Effect_Z_summary.txt", header=T, stringsAsFactor=F)
data_f<-data[grep("K562", data$CELL),]

data_f<-data.frame(data_f, P=pnorm(data_f$Z, lower.tail=FALSE, mean=0, sd=1))

traits<-read.table("RBC_traits_LDSC.txt", header=F,  sep="\t", stringsAsFactor=F)
row.names(traits)<-traits[,1]
data_f<-data_f[order(data_f$P),]
data_f<-data.frame(data_f, GROUP=ifelse(is.element(data_f$TRAIT, traits[,1]), "RBC traits", "others"), LABEL=traits[as.character(data_f$TRAIT), 2])
THRESH=0.05/nrow(data_f)

data_f$LABEL[data_f$P>THRESH]<-NA
data_f$TRAIT<-factor(as.character(data_f$TRAIT), levels=data_f$TRAIT)

library(ggplot2)
library(ggrepel)
g<-ggplot(data_f, aes(x=factor(TRAIT), y=-log10(P), color=GROUP, label=LABEL))
g<- g +  theme_classic()
g<- g + geom_point( shape=16, size=3 )
g<-g + theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), legend.title=element_blank())
g<-g + geom_hline(yintercept=-log10(THRESH), linetype="dashed")
g<-g + scale_x_discrete(expand=c(0.05,0))
g<-g + scale_color_manual(values = c("red", "black"), breaks=c("RBC traits", "others"))
g<-g + geom_text_repel(show.legend=FALSE, nudge_x = 10, direction = "y", hjust="left")
g<-g + ylab("-log10 (p value)")
g<-g + theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20), legend.text=element_text(size=20))
ggsave(plot=g, paste0("figure/Fig1B.pdf"), width=12, height=6)

##2. comparison of p values between K562 and MEP, RBC traits colored.

data<-read.table("result_summary/Effect_Z_summary.txt", header=T, stringsAsFactor=F)
data_f<-data[grep("K562|MEP", data$CELL),]
data_f<-data.frame(data_f, P=pnorm(data_f$Z, lower.tail=FALSE, mean=0, sd=1))

library(reshape2)
df<-dcast(data_f, TRAIT~CELL, value.var="P")

THRESH=0.05/nrow(df)

traits<-read.table("RBC_traits_LDSC.txt", header=F,  sep="\t", stringsAsFactor=F)
row.names(traits)<-traits[,1]
df<-df[order(df$MEP),]
df<-data.frame(df, GROUP=ifelse(is.element(df$TRAIT, traits[,1]), "RBC traits", "others"))

library(ggplot2)
library(ggrepel)
g<-ggplot(df, aes(x=-log10(MEP), y=-log10(K562), color=GROUP))
g<- g +  theme_classic()
g<- g + geom_point( shape=16, size=3 )
g<-g + theme(legend.title=element_blank())
g<-g + geom_hline(yintercept=-log10(THRESH), linetype="dashed")
g<-g + geom_vline(xintercept=-log10(THRESH), linetype="dashed")
g<-g + scale_color_manual(values = c("red", "black"), breaks=c("RBC traits", "others"))
g<-g + ylab("-log10 (p value), K562 enrichment") + xlab("-log10 (p value), MEP enrichment")
g<-g + theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20), legend.text=element_text(size=20))
ggsave(plot=g, paste0("figure/Fig1E.pdf"), width=8, height=6)


##3. heatmap of traits/cell types.
data<-read.table("result_summary/Effect_Z_summary.txt", header=T, stringsAsFactor=F)
data_f<-data
data_f<-data.frame(data_f, logP=-log10(pnorm(data_f$Z, lower.tail=FALSE, mean=0, sd=1)))

library(reshape2)
df<-dcast(data_f, CELL~TRAIT, value.var="logP")
row.names(df)<-df[,1]
df<-df[,-1]

granulo_traits<-c(30130, 30140, 30190, 30200, 30210, 30220 )
Platelet_traits<-c(30080,30090,30100,30110 )


library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

col_fun = colorRamp2(c(0,5,10 ), c("white", "orange", "red"))
ha = HeatmapAnnotation(
    RBCtraits = ifelse(is.element(colnames(df), traits[,1]), 1, 0),
    Mono_granulo_traits = ifelse(is.element(colnames(df), granulo_traits), 1, 0),
    Platelet_traits = ifelse(is.element(colnames(df), Platelet_traits), 1, 0),
    col = list(RBCtraits = c("1" = "blue", "0" = "white"), Mono_granulo_traits = c("1" = "blue", "0" = "white"),  Platelet_traits = c("1" = "blue", "0" = "white")),
    border = TRUE, show_legend=FALSE
)
pdf("figure/FigS1A.pdf", height=6, width=16)
set.seed(50)
ht=Heatmap(df, name="-log10 (p value)",
        row_names_gp = gpar(fontsize = 15), show_column_names=FALSE,
        col=col_fun,
        cluster_columns=TRUE, 
         show_column_dend = FALSE,
          cluster_rows=TRUE,show_row_dend = TRUE, row_names_side = "left",
        row_title=NULL, column_title=NULL,
         top_annotation=ha, 
 ) 
ht2=draw(ht, heatmap_legend_side = "left", annotation_legend_side = "right") 
ht2
dev.off()
