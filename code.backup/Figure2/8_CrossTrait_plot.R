summary<-data.frame()
for(TRAIT in c(86, 88)){
data<-read.table(paste0("data/LoF/GeneBayes_posterior/Backman_2021_", TRAIT, ".per_gene_estimates.tsv"), sep="\t", quote="", header=T, stringsAsFactor=F)
data<-data.frame(TRAIT=TRAIT, data)
summary<-rbind(summary, data)
}
summary$TRAIT<-gsub("86", "MCH", summary$TRAIT)
summary$TRAIT<-gsub("88", "RDW", summary$TRAIT)
summary$TRAIT<-gsub("106", "IRF", summary$TRAIT)

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

matt2<-t(mat2)
matt2<-data.frame(matt2, LABEL=row.names(matt2))

pca <- prcomp(cbind(matt2$MCH, matt2$RDW))$rotation
pca.slope <- pca[2,1] / pca[1,1]
pca.intercept <- mean(matt2$RDW) - (pca.slope * mean(matt2$MCH))

##color background 
df2=matt2
df2$LABEL[abs(df2$MCH)<0.45 & abs(df2$RDW)<0.45 ]<-""

tmp1=data.frame(ymin=min(df2$RDW)-0.1, ymax=0, xmin=min(df2$MCH), xmax=0, fill=1)
tmp2=data.frame(ymin=0, ymax=max(df2$RDW), xmin=0, xmax=max(df2$MCH), fill=1)
df.rect=rbind(tmp1, tmp2)
g<-ggplot(df2)
g<-g +  theme_classic(base_size = 40, base_family = "Helvetica")
g<-g + geom_rect(data=df.rect, aes(xmin=xmin, xmax=xmax, ymin=ymin,ymax=ymax), fill="mistyrose1")
g<-g + geom_point(aes(x=MCH, y=RDW), alpha=0.8, size=5)
g<-g + geom_hline(yintercept=0, linetype="dashed")
g<-g + geom_vline(xintercept=0, linetype="dashed")
g<-g + geom_text_repel(aes(x=MCH, y=RDW, label=LABEL), size=8,box.padding = 0.5, max.overlaps = Inf, point.size=3, min.segment.length=1 )
g<-g + geom_abline(aes(slope=pca.slope,intercept=pca.intercept),colour="black")
g<-g + xlab("MCH effect size (posterior)") + ylab("RDW effect size (posterior)")
ggsave(plot=g, "Fig2D.pdf", width=12,height=12)



