args <- commandArgs(trailingOnly = TRUE)
K <- if (length(args) >= 1) as.numeric(args[1]) else 60
result_dir <- if (length(args) >= 2) as.character(args[2]) else "program_trans/result"
regulator_file <- if (length(args) >= 3) as.character(args[3]) else paste0("data/Perturbseq/trait_association/K562GW/ProgramLevel/regulators_enrichment_K", K, "_Backman_2021_86.per_gene_estimates.tsv")
trait_label <- if (length(args) >= 4) as.character(args[4]) else "MCH"

regburden<-data.frame()

data<-read.table(regulator_file, header=T, stringsAsFactor=F)
data<-data.frame(Program=paste0("P", data$Program), data)
regburden<-rbind(regburden, data)

row.names(regburden)<-regburden[,1]
regburden<-data.frame(regburden, LOFscore=sign(regburden$beta_withShet)* (-log10(regburden$P_withShet)))

GWAS<-data.frame()
for(P in c(1:K)){
data<-read.table(file.path(result_dir, paste0("P", P, "_transeQTL_", trait_label, "_directionalTest.txt")), header=T, stringsAsFactor=F)
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
ggsave(plot=g, paste0(trait_label, "_transeQTLvsLoFburden_Fig4j.pdf"), width=12, height=12)
