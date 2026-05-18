
summary<-data.frame()
for(P in c(1:60)){
TF<-read.table(paste0("data/TF_ChIP/cNMF_enrichment/P", P, "_ChIPScoreenrichment.txt"), header=T, sep="\t", stringsAsFactor=F)
summary<-rbind(summary, TF)
}

summary<-data.frame(summary, FDR=p.adjust(summary$TFscore_program_MWU_P, method="BH"))

library(dplyr)
summary2<-data.frame()
for(P in c(1:60)){
tmp<-summary[is.element(summary$Program, P),]
tmp<- tmp %>% filter(FDR<0.05) %>% filter( TF_KD_Z < -2 & TF_type == "activator" |  TF_KD_Z > 2 & TF_type == "repressor" |  abs(TF_KD_Z) > 2 & TF_type == "both" )
summary2<-rbind(summary2, tmp)
}

##make heatmap
library(reshape2)
mat1<-dcast(summary, TF~Program, value.var="TFscore_program_MWU_P")
row.names(mat1)<-mat1[,1]
mat1<-mat1[,-1]
colnames(mat1)<-paste0("P", colnames(mat1))
mat1<--log10(mat1)
mat2<-matrix(rep(0, nrow(mat1)*ncol(mat1)), nrow=nrow(mat1))
row.names(mat2)<-row.names(mat1)
colnames(mat2)<-colnames(mat1)
for(i in 1:nrow(summary2)){
A<-summary2[i,"TF"]
B<-paste0("P", summary2[i,"Program"])
mat2[as.character(A), as.character(B)]<-1
}

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
col_fun = colorRamp2(c(0, -log10(0.05), -log10(0.02), 30 ), c("white", "white", "yellow", "red"))

mat3<-gsub(1, "\\*", mat2)
mat3<-gsub(0, "", mat3)
mat3<-matrix(mat3, nrow=nrow(mat2))
pdf("ChIP_score_enrich_program_FigS5.pdf", height=18, width=16)

set.seed(50)
ht=Heatmap(mat1, name="-lg10P",
        row_names_gp = gpar(fontsize = 10), show_column_names=TRUE, column_names_gp = gpar(fontsize = 12),
        col=col_fun,
        cluster_columns=T, 
         show_column_dend = T,
          cluster_rows=T,show_row_dend = T, row_names_side = "left",
        row_title=NULL, column_title=NULL,
         rect_gp = gpar(col = "black", lwd = 1),
         layer_fun =  function(j, i, x, y, w, h, fill) {
        v = pindex(mat3, i, j)
        gb = textGrob("*")
        gb_w = convertWidth(grobWidth(gb), "mm")
        gb_h = convertHeight(grobHeight(gb), "mm")
        l = v == "*"
        grid.text(sprintf("%s", v[l]), x[l], y[l] - gb_h*0.5 + gb_w*0.4, gp = gpar(fontsize = 10))
    }
 ) 
ht2=draw(ht, heatmap_legend_side = "right", annotation_legend_side = "left") 
ht2
dev.off()

