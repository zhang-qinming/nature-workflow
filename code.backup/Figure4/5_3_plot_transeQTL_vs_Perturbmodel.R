regburden<-data.frame()


data<-read.table("data/Perturbseq/trait_association/K562GW/ProgramLevel/regulators_enrichment_K60_Backman_2021_86.per_gene_estimates.tsv", header=T, stringsAsFactor=F)
data<-data.frame(Program=paste0("P", data$Program), data)
regburden<-rbind(regburden, data)

row.names(regburden)<-regburden[,1]
regburden<-data.frame(regburden, LOFscore=sign(regburden$beta_withShet)* (-log10(regburden$P_withShet)))


GWAS<-data.frame()
for(P in c(1:60)){
data<-read.table(paste0("result_mean/P", P, "_transeQTL_MCH_directionalTest.txt"), header=T, stringsAsFactor=F)
hoge<-data.frame(Program=paste0("P", P), GWASscore=sign(data$mean_Z_GWAS - data$mean_Z_ctrl)* (-log10(data$Ttest_P)))
GWAS<-rbind(GWAS, hoge)
}


df<-merge(regburden, GWAS, by="Program")
cor.test(df$LOFscore, df$GWASscore, method="p")
df<-data.frame(df, LABEL=df$Program)
df$LABEL[!is.element(df$LABEL, c("P16", "P40", "P4", "P25"))]<-""
df$LABEL<-gsub("P40", "Hemoglobin synthesis", df$LABEL)
df$LABEL<-gsub("P16", "Autophagy", df$LABEL)
df$LABEL<-gsub("P4", "S phase", df$LABEL)
df$LABEL<-gsub("P25", "G2/M phase", df$LABEL)

library(ggplot2)
library(ggrepel)

g<-ggplot(df, aes(x=LOFscore, y=GWASscore, label=LABEL))
g<-g +  theme_classic(base_size = 30, base_family = "Helvetica")
g<-g + geom_point(size=7)
g<-g + geom_vline(xintercept=0)
g<-g + geom_hline(yintercept=0)
g<-g + geom_text_repel(size=10)
g<-g + xlab("LOF regulator-burden correlation, signed -log10(p)")
g<-g + ylab("GWAS transeQTL effects, signed -log10(p)")
g<-g +theme(axis.text=element_text(size=35))
ggsave(plot=g, "MCH_transeQTLvsLoFburden_Fig4j.pdf", width=12, height=12)
