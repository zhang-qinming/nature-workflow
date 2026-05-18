
file_list=dir("data/Perturbseq/trait_association/K562GW/GeneLevel/")
file_list=file_list[grep("_geneRegulation_correlation.txt", file_list)]
file_list2<-file_list[grep("_86|_88|_96|_93|_101|_106", file_list)]

##
library(ggplot2)

summary3<-data.frame()
for(FILE in file_list2){
df<-read.table(paste0("data/BurdenRegCor/GeneLevel/", FILE), header=T, stringsAsFactor=F)
trait<-gsub(".per_gene_estimates.tsv_geneRegulation_correlation.txt", "", FILE)

df<-data.frame(df, signedP=sign(df$beta_withShet)*(-log10(df$P_withShet)))
df<-df[order(df$signedP),]
if(nrow(df) %% 2 ==0 ){
df<-data.frame(df, expected=sort(c( -log10(seq(1/nrow(df), 1, by=2/nrow(df))), log10(seq(1/nrow(df), 1, by=2/nrow(df))) )))
} else {
df<-data.frame(df, expected=sort(c( -log10(seq(1/nrow(df), 1, by=2/nrow(df))), log10(seq(2/nrow(df), 1, by=2/nrow(df))) )))	
}
colnames(df)<-gsub("P_withShet", "observed", colnames(df))
trait<-gsub("_mineto7", "", trait)
hoge<-data.frame(trait=trait, df)
summary3<-rbind(summary3, hoge)
}


library(ggplot2)
library(ggsci)
library(RColorBrewer)

corresp<-read.table("data/Perturbseq/trait_association/K562GW/GeneLevel/Backman_Niele_corresp_all", header=T, stringsAsFactor=F, sep="\t", quote="")
row.names(corresp)<-corresp$LOF_ID
summary3<-data.frame(summary3, TRAIT=corresp[ paste0(  summary3$trait, "_M1_001"), "Trait_GWAS"])
summary3$TRAIT[grep("_96", summary3$trait)] <- "Eosino count"
summary3$TRAIT<-gsub("Red blood cell \\(erythrocyte\\) distribution width", "RDW", summary3$TRAIT)
summary3$TRAIT<-gsub("Mean corpuscular haemoglobin", "MCH", summary3$TRAIT)
summary3$TRAIT<-gsub("Immature reticulocyte fraction", "IRF", summary3$TRAIT)
summary3$TRAIT<-gsub("Basophill percentage", "Basophil percentage", summary3$TRAIT)

summary3$TRAIT<-factor(as.character(summary3$TRAIT), levels=c("RDW", "MCH","IRF",  "Eosino count", "Basophil percentage",  "Lymphocyte count"))

##truncate the axis

summary4<-summary3[order(summary3$TRAIT, decreasing=T),]

summary4$signedP[summary4$signedP>10]<-10 ##truncate
summary4$signedP[summary4$signedP< -10]<- -10 ##truncate

g<-ggplot(summary4, aes(x=expected, y=signedP))
g<-g +  theme_set(theme_classic(base_family = "Helvetica", base_size=20))
g<-g + geom_point(aes(color=factor(TRAIT)), shape=16, size=3 )
g<-g + scale_color_manual(values = c(pal_igv("default")(15)[1:3], brewer.pal(n = 4, name = "Greys")[2:4]))
g<-g + geom_abline(intercept=0, slope=1, linetype="dashed")
g<-g + theme(legend.title=element_blank())
g<-g + xlab("E[Burden-regulator correlation, signed -log10(P)]") + ylab("O[Burden-regulator correlation, signed -log10(P)]")
ggsave(plot=g, paste0("Fig3D.pdf"), width=8, height=8)

