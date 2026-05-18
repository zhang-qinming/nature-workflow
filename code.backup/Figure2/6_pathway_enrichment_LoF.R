args   <- commandArgs(trailingOnly = T)   # enable reading factors from linux command line

trait <- as.character( args[1] ) ##86 for MCH
trait_GWAS<-as.character( args[1] ) ##30050 for MCH

traits<-read.table("data/Backman_Niele_corresp", header=T, stringsAsFactor=F)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)

LOF<-read.table(paste0("data/LoF/GeneBayes_posterior/Backman_2021_", trait, ".per_gene_estimates.tsv"), sep="\t", quote="", header=T, stringsAsFactor=F)
LOF<-LOF[order(abs(LOF$post_mean), decreasing=T),]

IDconv<-read.table("data/gencode_v41_gname_gid_ALL_sorted_onlyID", header=F, stringsAsFactor=F, row.names=1)
LOF<-data.frame(LOF, GENE=IDconv[LOF$ensg,1])
LOF<-na.omit(LOF)

GWAS<-read.table(paste0("data/GWAS/closest_genes/", trait_GWAS, "_genes_closest.txt"), header=F, stringsAsFactor=F)
N<-length(unique(GWAS[,1]))

ALLg<-LOF[c(1:N),"GENE"]

##GO
summary<-data.frame()


gse <- enrichGO(gene=ALLg, 
             ont ="ALL", 
             keyType = "SYMBOL",  
             minGSSize = 20, 
             maxGSSize = 2000, 
             qvalueCutoff = 1.5,
             pvalueCutoff = 1.5,
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "BH")

if(!is.null(gse)){
hoge<-data.frame(TOP_N="ALL", as.data.frame(gse))
summary<-rbind(summary, hoge)
}

trait_name=traits[is.element(traits$Backman, trait), "trait"]

write.table(data.frame(ID=row.names(summary), summary), paste0("enrichment_result/", trait_name, "_",  "LOF_GB_GO.txt"), row.names=F, sep="\t", quote=F)

##MsigDb
summary2<-data.frame()
library(msigdbr)

for(CAT in c("H")){

m_t2g <- msigdbr(species = "Homo sapiens", category = CAT) %>% 
  dplyr::select(gs_name, entrez_gene)

tmp2_entrez<-clusterProfiler:: bitr(ALLg, fromType = 'SYMBOL', toType = 'ENTREZID',OrgDb = org.Hs.eg.db, drop = TRUE)
tmp2_entrez<-tmp2_entrez[!duplicated(tmp2_entrez[,1]),]
tmp2_entrez<-unique(tmp2_entrez)[,2]

em <- enricher(tmp2_entrez, TERM2GENE=m_t2g, pvalueCutoff = 1.5,  qvalueCutoff=1.5)

if(nrow(em)>0){
hoge<-data.frame(TOP_N="ALL", as.data.frame(em))
summary2<-rbind(summary2, hoge)
}
}

write.table(data.frame(ID=row.names(summary2), summary2), paste0("enrichment_result/", trait_name, "_",  "LOF_GB_MsigDb.txt"), row.names=F, sep="\t", quote=F)
