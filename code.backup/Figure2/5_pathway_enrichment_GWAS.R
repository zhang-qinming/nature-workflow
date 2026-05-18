args   <- commandArgs(trailingOnly = T) 

trait <- as.character( args[1] ) ##loop over the traits
traits<-read.table("data/Backman_Niele_corresp", header=T, stringsAsFactor=F)

library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)
dir.create("enrichment_result", showWarnings=F)

GWAS<-read.table(paste0("data/GWAS/closest_genes/", trait, "_genes_closest.txt"), header=F, stringsAsFactor=F)
GWAS<-GWAS[order(GWAS[,2]),]
ALLg<-GWAS[,1]

##GO
summary<-data.frame()

gse <- enrichGO(gene=ALLg, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             minGSSize = 20, 
             maxGSSize = 2000, 
             pvalueCutoff = 1.5, 
             qvalueCutoff = 1.5, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "BH")


if(nrow(gse)>0){
hoge<-data.frame(TOP_N="ALL", as.data.frame(gse))
summary<-rbind(summary, hoge)
}


trait_name=traits[is.element(traits$Niele, trait), "trait"]

write.table(data.frame(ID=row.names(summary), summary), paste0("enrichment_result/", trait_name, "_GWAS_GO.txt"), row.names=F, sep="\t", quote=F)

##MsigDB
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

write.table(data.frame(ID=row.names(summary2), summary2), paste0("enrichment_result/", trait_name, "_GWAS_MsigDb.txt"), row.names=F, sep="\t", quote=F)
