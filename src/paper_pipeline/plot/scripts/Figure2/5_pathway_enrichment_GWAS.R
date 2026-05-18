args   <- commandArgs(trailingOnly = T) 

trait <- as.character( args[1] ) ##loop over the traits
traits<-read.table("data/Backman_Niele_corresp", header=T, stringsAsFactor=F)

library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)
script_path_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_path_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_path_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  getwd()
}
source(file.path(script_dir, "enrichment_helpers.R"))
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


if(has_rows(gse)){
hoge<-data.frame(TOP_N="ALL", as.data.frame(gse))
summary<-rbind(summary, hoge)
}


trait_name=traits[is.element(traits$Niele, trait), "trait"]

write_enrichment_table(summary, paste0("enrichment_result/", trait_name, "_GWAS_GO.txt"))

##MsigDB
summary2<-data.frame()
m_t2g <- load_hallmark_term2gene()

tmp2_entrez<-clusterProfiler:: bitr(ALLg, fromType = 'SYMBOL', toType = 'ENTREZID',OrgDb = org.Hs.eg.db, drop = TRUE)
tmp2_entrez<-tmp2_entrez[!duplicated(tmp2_entrez[,1]),]
tmp2_entrez<-unique(tmp2_entrez)[,2]

em <- enricher(tmp2_entrez, TERM2GENE=m_t2g, pvalueCutoff = 1.5,  qvalueCutoff=1.5)

if(has_rows(em)){
hoge<-data.frame(TOP_N="ALL", as.data.frame(em))
summary2<-rbind(summary2, hoge)
}

write_enrichment_table(summary2, paste0("enrichment_result/", trait_name, "_GWAS_MsigDb.txt"))
