
df<-read.table("multiCell_model/Regulators_correlation.txt",header=T, stringsAsFactor=F)

library(reshape2)
mat<-dcast(df, A~B, value.var="R")
row.names(mat)<-mat[,1]
mat<-mat[,-1]
mat[is.na(mat)]<-1

library(pheatmap)
pdf("regulators_correlation_heatmap.pdf")
p<-pheatmap(mat)
p
dev.off()

ORDER=row.names(mat)[p$tree_row$order]

##
burden_1<-data.frame()

CT_list=dir("data/Perturbseq/trait_association/")
for(CT in CT_list){
file_list=dir(paste0("data/Perturbseq/trait_association/", CT, "/ProgramLevel"))
file_list=file_list[grep("86", file_list)]
file_list=file_list[grep("regulators", file_list)]
FILE=file_list[grep("Backman", file_list)]

data<-read.table(paste0("data/Perturbseq/trait_association/", CT, "/ProgramLevel/", FILE), header=T, stringsAsFactor=F)
data<-data.frame(ID=paste0(CT, "_P", data$Program), data)
burden_1<-rbind(burden_1, data)
}
row.names(burden_1)<-burden_1[,1]
burden_1<-burden_1[colnames(mat),]

##
burden_2<-data.frame()

for(CT in CT_list){
file_list=dir(paste0("data/Perturbseq/trait_association/", CT, "/ProgramLevel"))
file_list=file_list[grep("88", file_list)]
file_list=file_list[grep("regulators", file_list)]
FILE=file_list[grep("Backman", file_list)]

data<-read.table(paste0("data/Perturbseq/trait_association/", CT, "/ProgramLevel/", FILE), header=T, stringsAsFactor=F)
data<-data.frame(ID=paste0(CT, "_P", data$Program), data)
burden_2<-rbind(burden_2, data)
}
row.names(burden_2)<-burden_2[,1]
burden_2<-burden_2[colnames(mat),]

##
burden_3<-data.frame()

for(CT in CT_list){
file_list=dir(paste0("data/Perturbseq/trait_association/", CT, "/ProgramLevel"))
file_list=file_list[grep("106", file_list)]
file_list=file_list[grep("regulators", file_list)]
FILE=file_list[grep("Backman", file_list)]

data<-read.table(paste0("data/Perturbseq/trait_association/", CT, "/ProgramLevel/", FILE), header=T, stringsAsFactor=F)
data<-data.frame(ID=paste0(CT, "_P", data$Program), data)
burden_3<-rbind(burden_3, data)
}
row.names(burden_3)<-burden_3[,1]
burden_3<-burden_3[colnames(mat),]

##
burden_4<-data.frame()

for(CT in CT_list){
file_list=dir(paste0("data/Perturbseq/trait_association/", CT, "/ProgramLevel"))
file_list=file_list[grep("131", file_list)]
file_list=file_list[grep("regulators", file_list)]
FILE=file_list[grep("Backman", file_list)]

data<-read.table(paste0("data/Perturbseq/trait_association/", CT, "/ProgramLevel/", FILE), header=T, stringsAsFactor=F)
data<-data.frame(ID=paste0(CT, "_P", data$Program), data)
burden_4<-rbind(burden_4, data)
}

row.names(burden_4)<-burden_4[,1]
burden_4<-burden_4[colnames(mat),]

##
burden_5<-data.frame()

for(CT in CT_list){
file_list=dir(paste0("data/Perturbseq/trait_association/", CT, "/ProgramLevel"))
file_list=file_list[grep("132", file_list)]
file_list=file_list[grep("regulators", file_list)]
FILE=file_list[grep("Backman", file_list)]

data<-read.table(paste0("data/Perturbseq/trait_association/", CT, "/ProgramLevel/", FILE), header=T, stringsAsFactor=F)
data<-data.frame(ID=paste0(CT, "_P", data$Program), data)
burden_5<-rbind(burden_5, data)
}

row.names(burden_5)<-burden_5[,1]
burden_5<-burden_5[colnames(mat),]

##
burden_6<-data.frame()

for(CT in CT_list){
file_list=dir(paste0("data/Perturbseq/trait_association/", CT, "/ProgramLevel"))
file_list=file_list[grep("93", file_list)]
file_list=file_list[grep("regulators", file_list)]
FILE=file_list[grep("Backman", file_list)]

data<-read.table(paste0("data/Perturbseq/trait_association/", CT, "/ProgramLevel/", FILE), header=T, stringsAsFactor=F)
data<-data.frame(ID=paste0(CT, "_P", data$Program), data)
burden_6<-rbind(burden_6, data)
}

row.names(burden_6)<-burden_6[,1]
burden_6<-burden_6[colnames(mat),]

##
burden_7<-data.frame()

for(CT in CT_list){
file_list=dir(paste0("data/Perturbseq/trait_association/", CT, "/ProgramLevel"))
file_list=file_list[grep("133", file_list)]
file_list=file_list[grep("regulators", file_list)]
FILE=file_list[grep("Backman", file_list)]

data<-read.table(paste0("data/Perturbseq/trait_association/", CT, "/ProgramLevel/", FILE), header=T, stringsAsFactor=F)
data<-data.frame(ID=paste0(CT, "_P", data$Program), data)
burden_7<-rbind(burden_7, data)
}

row.names(burden_7)<-burden_7[,1]
burden_7<-burden_7[colnames(mat),]


##
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

col_fun = colorRamp2(c(-1, 0, 1 ), c("blue", "white", "red"))

tmp<-rbind(burden_1, burden_2, burden_3,burden_4,burden_5, burden_6,burden_7)

tmp2<-data.frame(P=tmp$P_withShet, FDR=p.adjust(tmp$P_withShet, method="BH"))
tmp2<-tmp2[order(tmp2$P),]
THRESH=tmp2$P[tmp2$FDR>0.01][1] ##set FDR to 0.01

burden_1<-data.frame(burden_1, COLOR=ifelse( (burden_1$P_withShet<THRESH & burden_1$beta_withShet >0), "red", ifelse ( (burden_1$P_withShet<THRESH & burden_1$beta_withShet <0), "blue", "grey")))
burden_2<-data.frame(burden_2, COLOR=ifelse( (burden_2$P_withShet<THRESH & burden_2$beta_withShet >0), "red", ifelse ( (burden_2$P_withShet<THRESH & burden_2$beta_withShet <0), "blue", "grey")))
burden_3<-data.frame(burden_3, COLOR=ifelse( (burden_3$P_withShet<THRESH & burden_3$beta_withShet >0), "red", ifelse ( (burden_3$P_withShet<THRESH & burden_3$beta_withShet <0), "blue", "grey")))
burden_4<-data.frame(burden_4, COLOR=ifelse( (burden_4$P_withShet<THRESH & burden_4$beta_withShet >0), "red", ifelse ( (burden_4$P_withShet<THRESH & burden_4$beta_withShet <0), "blue", "grey")))
burden_5<-data.frame(burden_5, COLOR=ifelse( (burden_5$P_withShet<THRESH & burden_5$beta_withShet >0), "red", ifelse ( (burden_5$P_withShet<THRESH & burden_5$beta_withShet <0), "blue", "grey")))
burden_6<-data.frame(burden_6, COLOR=ifelse( (burden_6$P_withShet<THRESH & burden_6$beta_withShet >0), "red", ifelse ( (burden_6$P_withShet<THRESH & burden_6$beta_withShet <0), "blue", "grey")))
burden_7<-data.frame(burden_7, COLOR=ifelse( (burden_7$P_withShet<THRESH & burden_7$beta_withShet >0), "red", ifelse ( (burden_7$P_withShet<THRESH & burden_7$beta_withShet <0), "blue", "grey")))

cell_type=ifelse(grepl("K562GW", colnames(mat)), "K562GW", ifelse(grepl("K562_essential", colnames(mat)), "K562essential", ifelse(grepl("hepg2", colnames(mat)), "HepG2", ifelse(grepl("jurkat", colnames(mat)), "Jurkat", 
    ifelse(grepl("GSM6858449", colnames(mat)),"THP1", "RPE1")))))


column_ha3 = HeatmapAnnotation(
  "MCH" = anno_barplot(-log10(burden_1$P_withShet),  gp = gpar(fill = burden_1$COLOR), axis_param=list(gp=gpar(fontsize = 12))), 
  "RDW" = anno_barplot(-log10(burden_2$P_withShet),  gp = gpar(fill = burden_2$COLOR), axis_param=list(gp=gpar(fontsize = 12))), 
  "IRF" = anno_barplot(-log10(burden_3$P_withShet),  gp = gpar(fill = burden_3$COLOR), axis_param=list(gp=gpar(fontsize = 12))), 
  "HBA1c" = anno_barplot(-log10(burden_4$P_withShet),  gp = gpar(fill = burden_4$COLOR), axis_param=list(gp=gpar(fontsize = 12))), 
  "HDL" = anno_barplot(-log10(burden_5$P_withShet),  gp = gpar(fill = burden_5$COLOR), axis_param=list(gp=gpar(fontsize = 12))), 
  "Lymphocyte\ncount" = anno_barplot(-log10(burden_6$P_withShet),  gp = gpar(fill = burden_6$COLOR), axis_param=list(gp=gpar(fontsize = 12))), 
  "IGF-1" = anno_barplot(-log10(burden_7$P_withShet),  gp = gpar(fill = burden_7$COLOR), axis_param=list(gp=gpar(fontsize = 12))), 
  "K562GW" = ifelse(cell_type=="K562GW", 1, 0), "K562essential" = ifelse(cell_type=="K562essential", 1, 0),
   "HepG2" = ifelse(cell_type=="HepG2", 1, 0),
    "Jurkat" = ifelse(cell_type=="Jurkat", 1, 0),
     "RPE1" = ifelse(cell_type=="RPE1", 1, 0),
  col = list("K562GW" = c("1" = "black", "0"="white"), "K562essential" = c("1" = "black", "0"="white"), "HepG2" = c("1" = "black", "0"="white"), "Jurkat" = c("1" = "black", "0"="white"), "RPE1" = c("1" = "black", "0"="white")),
height =  unit(160, "mm"), show_legend = FALSE, annotation_name_gp= gpar(fontsize = 18) )

set.seed(50)
ht1=Heatmap(mat, name="Regulator correlation",rect_gp = gpar(type = "none"), 
        #row_names_gp = gpar(fontsize = 3), 
        show_column_names=FALSE, show_row_names=FALSE,
        #column_names_gp = gpar(fontsize = 3),
        col=col_fun,
        #cluster_columns=T, 
        column_order = ORDER, row_order=ORDER,
        # show_column_dend = T,
        #  cluster_rows=T,show_row_dend = T, 
        #row_names_side = "left",
        row_title=NULL, column_title="Gene expression programs", column_title_side = "bottom",
         #rect_gp = gpar(col = "black", lwd = 1),　
         cell_fun = function(j, i, x, y, w, h, fill) {
        if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6 ){
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
        }
        },
        heatmap_legend_param = list(
        title = "Regulator\nPearson's r", 
        labels_gp=gpar(fontsize = 15), title_gp=gpar(fontsize = 15)
        ),
         top_annotation = column_ha3, use_raster = TRUE, raster_quality = 1)
##

##
df2<-read.table("multiCell_model/Program_genes_Jaccard.txt",header=T, stringsAsFactor=F)

library(reshape2)
mat2<-dcast(df2, A~B, value.var="JI")
row.names(mat2)<-mat2[,1]
mat2<-mat2[,-1]
#mat[is.na(mat)]<-1
col_fun2 = colorRamp2(c(0, 0.2, 0.5 ), c("white", "orange", "red"))
mat2<-mat2[,colnames(mat)]
mat2<-mat2[row.names(mat),]

ht2=Heatmap(mat2, name="Jaccard index", rect_gp = gpar(type = "none"), 
        #row_names_gp = gpar(fontsize = 3), 
        show_column_names=FALSE, show_row_names=FALSE,
        #column_names_gp = gpar(fontsize = 3),
        col=col_fun2,
        #cluster_columns=T, 
        column_order = ORDER, row_order=ORDER,
        # show_column_dend = T,
        #  cluster_rows=T,show_row_dend = T, 
        #row_names_side = "left",
        row_title=NULL, column_title="Gene expression programs", column_title_side = "bottom",
          cell_fun = function(j, i, x, y, w, h, fill) {
        if(as.numeric(x) >= 1 - as.numeric(y) - 1e-6){
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
        }},
        heatmap_legend_param = list(
        title = "Program\nJaccard index", 
        labels_gp=gpar(fontsize = 15), title_gp=gpar(fontsize = 15)
        ), use_raster = TRUE, raster_quality = 1
        
         #rect_gp = gpar(col = "black", lwd = 1),　
         #top_annotation = column_ha3
         )



pdf("heatmap_withBurdenCor_all_programs_suppFigure.pdf", height=13, width=15)
set.seed(50)
draw(ht1 + ht2, ht_gap = unit(-336, "mm"))
dev.off()
