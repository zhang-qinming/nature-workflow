args   <- commandArgs(trailingOnly = T)  
options("scipen"=10) ##disable 100000 to 1e6
FILE<-as.character( args[1] ) ##e.g., Backman_2021_86.per_gene_estimates.tsv
K<-as.numeric( args[2] ) ##60

pro<-read.table(paste0("data/Perturbseq/trait_association/K562GW/ProgramLevel/programs_enrichment_K", K, "_", FILE), header=T, stringsAsFactor=F)
reg<-read.table(paste0("data/Perturbseq/trait_association/K562GW/ProgramLevel/regulators_enrichment_K", K, "_", FILE), header=T, stringsAsFactor=F)

df<-merge(pro, reg, by="Program")

df<-data.frame(df, Program_score= sign( df$MEANgamma_top100 - df$shet_adjusted_random_mean) * (-log10(df$MEANgamma_top100_shet_adjusted_P)), 
	regulator_score=sign(df$beta_withShet) * (-log10(df$P_withShet)))
df<-data.frame(df, LABEL=df$Program)
df$LABEL[df$MEANgamma_top100_shet_adjusted_P > 0.05/nrow(df) & df$P_withShet > 0.05/nrow(df) & !is.element(df$Program, stepwise_selected) ]<-""
df<-data.frame(df, COLOR=ifelse(df$MEANgamma_top100_shet_adjusted_P < 0.05/nrow(df) & df$P_withShet < 0.05/nrow(df), "both_enriched", 
	ifelse(df$MEANgamma_top100_shet_adjusted_P < 0.05/nrow(df) & df$P_withShet >= 0.05/nrow(df), "program_enriched",  
	ifelse(df$MEANgamma_top100_shet_adjusted_P >= 0.05/nrow(df) & df$P_withShet < 0.05/nrow(df), "regulator_enriched",  "other"))))

df$COLOR<-factor(as.character(df$COLOR), levels=c("other", "program_enriched", "regulator_enriched", "both_enriched"))
library(ggplot2)
library(ggrepel)
g<-ggplot(df, aes(x=Program_score, y=regulator_score, label=LABEL, color=COLOR))
g<-g +  theme_classic(base_size = 30, base_family = "Helvetica")
g<-g + geom_point(size=6)
g<-g + geom_vline(xintercept=0, linetype="dashed")
g<-g + geom_hline(yintercept=0, linetype="dashed")
g<-g + geom_text_repel(size=10)
g<-g + scale_color_manual(values=c("other"="grey", "program_enriched"="#FEA601", "regulator_enriched"="#4783B5", "both_enriched"="#34a853"))
g<-g + guides(color = guide_legend(override.aes = aes(label = "", alpha = 1)))
g<-g + theme(legend.title=element_blank())
g<-g + xlab("Program burden effect, signed -log10(P)") + ylab("Regulator-burden correlation, signed -log10(P)")
ggsave(plot=g, paste0(FILE, "_Fig4C.pdf"), width=15, height=10)
