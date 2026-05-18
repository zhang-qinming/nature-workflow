
args <- commandArgs(trailingOnly = TRUE)

corresp<-read.table("data/gencode_v41_gname_gid_ALL_sorted_onlyID", header=F, stringsAsFactor=F)
corresp<-corresp[!duplicated(corresp[,1]),]
row.names(corresp)<-corresp[,1]


##
trait <- if (length(args) >= 1) as.character(args[1]) else "30050"
COMP <- if (length(args) >= 2) as.character(args[2]) else "closest"
target_gene <- if (length(args) >= 3) as.character(args[3]) else "HBA1"
lof_burden_file <- if (length(args) >= 4) as.character(args[4]) else "Backman_2021_86_M1_001.summary_statistics.csv"
GWAS<-read.table(paste0("data/GWAS/closest_genes/", trait, "_genes_", COMP, ".txt"), header=F, stringsAsFactor=F)

LOF_rawZ<-read.csv(paste0("data/LoF/raw_burden/", lof_burden_file), header=T, stringsAsFactor=F)
LOF_rawZ<-data.frame(LOF_rawZ, Z=LOF_rawZ$beta/LOF_rawZ$standard_error)
LOF_rawZ<-LOF_rawZ[order(abs(LOF_rawZ$Z), decreasing=T),]
LOF_rawZ<-data.frame(GENE=corresp[LOF_rawZ$ensg, 2], LOF_rawZ)
LOF_rawZ<-LOF_rawZ[,c("GENE", "Z")]
LOF_rawZ<-na.omit(LOF_rawZ)

colnames(GWAS)<-c("GENE", "P")

##HBA1 regulators

library(data.table)
df1<-fread("data/Perturbseq/geneLevel/K562GW/limma_P_sum.txt", data.table=F, header=T)
row.names(df1)<-df1[,1]
df1<-df1[,-1]

corresp<-read.table("data/gencode_v41_gname_gid_ALL_sorted_onlyID", header=F, stringsAsFactor=F)
corresp<-corresp[!duplicated(corresp[,1]),]
corresp<-corresp[!duplicated(corresp[,2]),]
row.names(corresp)<-corresp[,2]
library("boot")

ensg<-corresp[is.element(corresp[,2], target_gene), 1]

df1_f<-df1[ensg,,drop=F]

##GWAS enrichment
gtf<-read.table("data/GWAS/genes.protein_coding.v39.gtf", header=T, stringsAsFactor=F)
colnames(gtf)[7]<-"GENE"

df1_f2<-t(df1_f)
df1_f2<-df1_f2[order(df1_f2[,1]),,drop=F]
df1_f2<-data.frame(df1_f2, FDR=p.adjust(df1_f2[,1], method="BH"))
BG_HBA1<-intersect(names(df1_f), gtf$GENE)

summary<-data.frame()

for(HBA1_N in c(200, 416, 605, 800, 1000)){
HBA1_regulators=row.names(df1_f2)[1:HBA1_N]
for(N in c(90, 543)){
topgenes<-GWAS$GENE[1:N]
A<-length(unique(intersect(HBA1_regulators, topgenes)))
tmp<-length(unique(intersect(HBA1_regulators, BG_HBA1)))
B<-tmp - A
tmp2<-length(unique(intersect(topgenes, BG_HBA1)))
C<-tmp2 - A
D<- length(BG_HBA1) - A - B - C

mat<-matrix(c(A,B,C,D), nrow=2, byrow=T)

fisher_P<-fisher.test(mat, alternative="two.sided")$p.value
fisher_OR<-fisher.test(mat, alternative="two.sided")$estimate
fisher_OR_lowCI<-fisher.test(mat, alternative="two.sided")$conf.int[1]
fisher_OR_highCI<-fisher.test(mat, alternative="two.sided")$conf.int[2]

hoge<-data.frame(TEST="GWAS", HBA1regN=HBA1_N, topN=N, COMP="HBA1regulators", fisher_P=fisher_P, fisher_OR=fisher_OR, fisher_OR_lowCI=fisher_OR_lowCI, fisher_OR_highCI=fisher_OR_highCI)
summary<-rbind(summary, hoge)
}
}

write.table(summary, paste0("MCH_GWAS_", target_gene, "regulator_enrichment.txt"), row.names=F, sep="\t", quote=F)


##LOF enrichment

BG_HBA1<-intersect(names(df1_f), LOF_rawZ$GENE)
summary<-data.frame()
for(HBA1_N in c(200, 416, 605, 800, 1000)){
HBA1_regulators=row.names(df1_f2)[1:HBA1_N]

for(N in c(90)){
topgenes<-LOF_rawZ$GENE[1:N]

##HBA1regulator enrichment
A<-length(unique(intersect(HBA1_regulators, topgenes)))
tmp<-length(unique(intersect(HBA1_regulators, BG_HBA1)))
B<-tmp - A
tmp2<-length(unique(intersect(topgenes, BG_HBA1)))
C<-tmp2 - A
D<- length(BG_HBA1) - A - B - C

mat<-matrix(c(A,B,C,D), nrow=2, byrow=T)

fisher_P<-fisher.test(mat, alternative="two.sided")$p.value
fisher_OR<-fisher.test(mat, alternative="two.sided")$estimate
fisher_OR_lowCI<-fisher.test(mat, alternative="two.sided")$conf.int[1]
fisher_OR_highCI<-fisher.test(mat, alternative="two.sided")$conf.int[2]

hoge<-data.frame(TEST="LOF",  HBA1regN=HBA1_N,  topN=N, COMP="HBA1regulators", fisher_P=fisher_P, fisher_OR=fisher_OR, fisher_OR_lowCI=fisher_OR_lowCI, fisher_OR_highCI=fisher_OR_highCI)
summary<-rbind(summary, hoge)
}
}

write.table(summary, paste0("MCH_LoF_", target_gene, "regulator_enrichment.txt"), row.names=F, sep="\t", quote=F)

##make figure

LOF_enrich<-read.table(paste0("MCH_LoF_", target_gene, "regulator_enrichment.txt"), header=T, stringsAsFactor=F)
GWAS_enrich<-read.table(paste0("MCH_GWAS_", target_gene, "regulator_enrichment.txt"), header=T, stringsAsFactor=F)

df<-rbind(LOF_enrich, GWAS_enrich)
df<-data.frame(df, TEST2=paste0(df$TEST, "_", df$topN))

library(ggplot2)


df$fisher_OR_highCI[df$fisher_OR_highCI > 10] <- 10  ##truncate
df$TEST2<-gsub("_", " top ", df$TEST2)
df$TEST2<-factor(as.character(df$TEST2), levels=c("LOF top 90", "GWAS top 90", "GWAS top 543"))

g1<-ggplot(df, aes(x=factor(HBA1regN), y=fisher_OR, color=TEST2))
g1<-g1 +  theme_classic(base_size = 30, base_family = "Helvetica")
g1<-g1 + geom_point(size=8, position = position_dodge(.6))
g1<-g1 + geom_errorbar(aes(ymax=fisher_OR_highCI, ymin=fisher_OR_lowCI), width=0,position=position_dodge(.6))
g1<-g1 + geom_hline(yintercept=1, linetype="dashed")
g1<-g1 + ylim(c(0,10))
g1<-g1 + ylab("odds ratio")
g1<-g1 + xlab("Top n of genes, HBA1 regulators")
g1<- g1 + scale_color_manual(values=c("#56941E", "#A874D2", "#652F90")) ##purple green
ggsave(plot=g1, "Fig3B.pdf", width=12, height=8)
