
args <- commandArgs(trailingOnly = TRUE)
trait <- if (length(args) >= 1) as.character(args[1]) else "30050"

library(data.table)

gwas<-fread(paste0("data/GWAS/", trait, "_irnt.tsv.gz"))
gwas_f<-gwas[gwas$low_confidence_variant == FALSE,]
gwas_f2<-data.frame(CHROM=sapply(strsplit(as.character(gwas_f$variant), ":"), function(x){x[1]}), 
	POS=as.numeric(sapply(strsplit(as.character(gwas_f$variant), ":"), function(x){x[2]})),
	P=gwas_f$pval)

gwas_f3<-data.frame(gwas_f2[is.element(gwas_f2$CHROM, c(1:22)),], COORD=0)
for(chr in c(1:22)){
if(chr==1){
gwas_f3$COORD[is.element(gwas_f3$CHROM, chr)]<- gwas_f3$POS[is.element(gwas_f3$CHROM, chr)]
} else {
gwas_f3$COORD[is.element(gwas_f3$CHROM, chr)]<- gwas_f3$POS[is.element(gwas_f3$CHROM, chr)] + max(  na.omit(gwas_f3$COORD[is.element(gwas_f3$CHROM, ( chr - 1 ))] ))
}
}


gtf<-read.table("data/GWAS/genes.protein_coding.v39.gtf", header=T, stringsAsFactor=F)

gwas_heme<-data.frame()
for(geneset in c("HALLMARK_HEME_METABOLISM", "Hematopoiesisgenes", "mitotic_cell_cycle", "positive_macromolecule_synthesis" )){
heme<-read.table(paste0("data/geneset/", geneset, ".txt"), header=F, stringsAsFactor=F)[,1]
heme_pos<-gtf[is.element(gtf$gene, heme), c("gene", "chr", "start", "end")]
heme_pos<-data.frame(heme_pos, P1=heme_pos$start - 50000, P2=heme_pos$end + 50000)

library(dplyr)
for(i in 1:nrow(heme_pos)){
hoge<- gwas_f3 %>% filter(CHROM==gsub("chr", "", heme_pos[i,2])) %>%   filter(POS > heme_pos[i,"P1"]) %>%   filter(POS < heme_pos[i,"P2"])  %>%   filter(P<5E-8)
if(nrow(hoge)>0){
hoge<-data.frame(hoge, LABEL="", geneset=geneset)
hoge<-hoge[order(hoge$P),]
hoge$LABEL[1]<-heme_pos[i,"gene"]
hoge<-hoge[order(hoge$POS),]
gwas_heme<-rbind(gwas_heme, hoge)
}
}
}
gwas_heme<-gwas_heme[!duplicated(gwas_heme$COORD),] ##remove duplicated rows. If >2 geneset were annotated, represent by the first gene set.

gwas_heme$LABEL[grep("^HB", gwas_heme$LABEL)] <- "Hemoglobin genes"
gwas_heme$LABEL[grep("^H4", gwas_heme$LABEL)] <- "H4 family genes"
gwas_heme$LABEL[grep("HLA", gwas_heme$LABEL)] <- "HLA genes"

gwas_heme<-gwas_heme[order(gwas_heme$P),]
gwas_heme$LABEL[duplicated(gwas_heme$LABEL)]<-""
gwas_heme$LABEL[gwas_heme$P>1E-30] <- ""
gwas_heme$P[gwas_heme$P==0] <- 1E-321

library(topr)
library(ggplot2)
library(ggrastr)
library(ggrepel)

gwas_f3<-gwas_f3[gwas_f3$P>0,]

library(dplyr)
axis_set <- gwas_f3 |>
  group_by(CHROM) |>
  summarize(center = mean(COORD))

library(scattermore)
library(ggsci)

gwas_f3<-data.frame(gwas_f3, geneset="other")


gwas_f3$geneset<-gsub("HALLMARK_HEME_METABOLISM", "Heme metabolism", gwas_f3$geneset)
gwas_f3$geneset<-gsub("Hematopoiesisgenes", "Hematopoiesis", gwas_f3$geneset)
gwas_f3$geneset<-gsub("mitotic_cell_cycle", "Mitotic cell cycle", gwas_f3$geneset)
gwas_f3$geneset<-gsub("positive_macromolecule_synthesis", "Macromolecule synthesis, positive reg", gwas_f3$geneset)

gwas_heme$geneset<-gsub("HALLMARK_HEME_METABOLISM", "Heme metabolism", gwas_heme$geneset)
gwas_heme$geneset<-gsub("Hematopoiesisgenes", "Hematopoiesis", gwas_heme$geneset)
gwas_heme$geneset<-gsub("mitotic_cell_cycle", "Mitotic cell cycle", gwas_heme$geneset)
gwas_heme$geneset<-gsub("positive_macromolecule_synthesis", "Macromolecule synthesis, positive reg", gwas_heme$geneset)


gwas_f4_1<-gwas_f3[gwas_f3$P<1E-3,]
gwas_f4_2<-gwas_f3[gwas_f3$P >= 1E-3,]
gwas_f4<-rbind(gwas_f4_1, gwas_f4_2)

g<-ggplot(gwas_f4, aes(x=COORD, y=-log10(P)))
g<-g +  theme_classic(base_size = 30, base_family = "Helvetica")
g<-g + geom_point_rast(alpha=0.5,  size=2, aes(color=geneset), shape=16)
g<-g + geom_hline(yintercept=-log10(5E-8), linetype="dashed", color="red")
g<-g + geom_point(data=gwas_heme, aes(color=geneset), size=2, shape=16)
g<-g + geom_text_repel(data=gwas_heme, aes(label=LABEL, color=geneset), size=8, nudge_y=10, max.overlaps=100, show.legend = FALSE )
g<-g + scale_color_manual(values = c("other"="grey60", "Heme metabolism"=pal_d3("category10")(4)[1], "Hematopoiesis"=pal_d3("category10")(4)[2],
  "Mitotic cell cycle"=pal_d3("category10")(4)[3], "Macromolecule synthesis, positive reg"=pal_d3("category10")(4)[4]),
breaks=c("Heme metabolism", "Hematopoiesis", "Mitotic cell cycle", "Macromolecule synthesis, positive reg", "other" ))
g<-g + scale_x_continuous(
    label = axis_set$CHROM,
    breaks = axis_set$center) 
g<-g + scale_y_continuous(expand = c(0, 0), limits=c(0,390))
g<-g + xlab("chromosome")
ggsave(plot=g, "Fig2A.png", width=30, height=12, dpi=600)
ggsave(plot=g, "Fig2A.pdf", width=30, height=12, dpi=600) ##partly rasterized. can be treated with illustrator
