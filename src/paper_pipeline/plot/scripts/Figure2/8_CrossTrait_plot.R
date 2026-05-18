args <- commandArgs(trailingOnly = TRUE)
trait_ids <- if (length(args) >= 1) strsplit(as.character(args[1]), ",")[[1]] else c("86", "88")
trait_ids <- trimws(trait_ids)
if (length(trait_ids) != 2) {
  stop("8_CrossTrait_plot.R expects exactly two comma-separated trait IDs")
}

trait_label <- function(trait_id) {
  if (trait_id == "86") {
    return("MCH")
  }
  if (trait_id == "88") {
    return("RDW")
  }
  if (trait_id == "106") {
    return("IRF")
  }
  paste0("Trait", trait_id)
}

trait_x <- trait_label(trait_ids[1])
trait_y <- trait_label(trait_ids[2])

summary<-data.frame()
for(TRAIT in trait_ids){
data<-read.table(paste0("data/LoF/GeneBayes_posterior/Backman_2021_", TRAIT, ".per_gene_estimates.tsv"), sep="\t", quote="", header=T, stringsAsFactor=F)
data<-data.frame(TRAIT=trait_label(TRAIT), data)
summary<-rbind(summary, data)
}

library(reshape2)
mat2<-dcast(summary, TRAIT~ensg, value.var="post_mean")
row.names(mat2)<-mat2[,1]
mat2<-mat2[,-1]

corresp<-read.table("data/gencode_v41_gname_gid_ALL_sorted_onlyID", header=F, stringsAsFactor=F)
corresp<-corresp[!duplicated(corresp[,1]),]
corresp<-corresp[!duplicated(corresp[,2]),]
row.names(corresp)<-corresp[,1]

colnames(mat2)<-corresp[colnames(mat2),2]


library(ggplot2)
library(ggrepel)

matt2<-t(mat2)
matt2<-data.frame(matt2, LABEL=row.names(matt2))

pca <- prcomp(cbind(matt2[[trait_x]], matt2[[trait_y]]))$rotation
pca.slope <- pca[2,1] / pca[1,1]
pca.intercept <- mean(matt2[[trait_y]]) - (pca.slope * mean(matt2[[trait_x]]))

##color background 
df2=matt2
df2$LABEL[abs(df2[[trait_x]])<0.45 & abs(df2[[trait_y]])<0.45 ]<-""

tmp1=data.frame(ymin=min(df2[[trait_y]])-0.1, ymax=0, xmin=min(df2[[trait_x]]), xmax=0, fill=1)
tmp2=data.frame(ymin=0, ymax=max(df2[[trait_y]]), xmin=0, xmax=max(df2[[trait_x]]), fill=1)
df.rect=rbind(tmp1, tmp2)
g<-ggplot(df2)
g<-g +  theme_classic(base_size = 40, base_family = "Helvetica")
g<-g + geom_rect(data=df.rect, aes(xmin=xmin, xmax=xmax, ymin=ymin,ymax=ymax), fill="mistyrose1")
g<-g + geom_point(aes_string(x=trait_x, y=trait_y), alpha=0.8, size=5)
g<-g + geom_hline(yintercept=0, linetype="dashed")
g<-g + geom_vline(xintercept=0, linetype="dashed")
g<-g + geom_text_repel(aes_string(x=trait_x, y=trait_y, label="LABEL"), size=8,box.padding = 0.5, max.overlaps = Inf, point.size=3, min.segment.length=1 )
g<-g + geom_abline(aes(slope=pca.slope,intercept=pca.intercept),colour="black")
g<-g + xlab(paste0(trait_x, " effect size (posterior)")) + ylab(paste0(trait_y, " effect size (posterior)"))
ggsave(plot=g, "Fig2D.pdf", width=12,height=12)
