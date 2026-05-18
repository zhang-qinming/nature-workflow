
args   <- commandArgs(trailingOnly = T)  
options("scipen"=10) ##disable 100000 to 1e6
dir.create("multiCell_model/conditional", showWarnings=F)

K<-60 ##60
cell_list=dir("data/Perturbseq/cNMF_regulation/")

GEP_reg_beta<-data.frame()

for(CELL in cell_list){
##data loading
for(i in 1:K){
tmp<-read.table(paste0("data/Perturbseq/cNMF_regulation/", CELL, "/K60_program", i, "_perturb_effects.txt"), header=T)
tmp1<-tmp[,c(1,2)]

colnames(tmp1)<-c("GENE", paste0(CELL, "_", "P", i))

if(i==1 & CELL==cell_list[1]){
GEP_reg_beta<-tmp1
} else {
GEP_reg_beta<-merge(GEP_reg_beta, tmp1, by="GENE")
}
}
}

##

FILE="Backman_2021_131.per_gene_estimates.tsv"

LOF<-read.table(paste0("data/LoF/GeneBayes_posterior/", FILE), sep="\t", quote="", header=T, stringsAsFactor=F)

LOF$post_mean[is.element(LOF$post_mean, "Inf")]<-max(LOF$post_mean[!is.infinite(LOF$post_mean)])
LOF$post_mean[is.element(LOF$post_mean, "-Inf")]<-min(LOF$post_mean[!is.infinite(LOF$post_mean)])

corresp<-read.table("data/gencode_v41_gname_gid_ALL_sorted_onlyID", header=F, stringsAsFactor=F)
corresp<-corresp[!duplicated(corresp[,1]),]
ensg_to_gene <- corresp[,2]
names(ensg_to_gene) <- corresp[,1]

LOF<-data.frame(LOF, gene=unname(ensg_to_gene[as.character(LOF$ensg)]))

shet<-read.table("data/shet_10bins.txt", header=T, stringsAsFactor=F)
shet<-shet[is.element(shet$ensg, LOF$ensg),]

corresp<-read.table("data/gencode_v41_gname_gid_ALL_sorted_onlyID", header=F, stringsAsFactor=F)
corresp<-corresp[!duplicated(corresp[,1]),]
corresp<-corresp[!duplicated(corresp[,2]),]
symbol_to_ensg <- corresp[,1]
names(symbol_to_ensg) <- corresp[,2]


summary_reg<-data.frame()

program_list=colnames(GEP_reg_beta)
program_list=program_list[-1]

top_list=c("GSE264667_jurkat_raw_singlecell_01_P42", "K562GW_P28", "rpe1_raw_singlecell_01_P21")
for(Program in program_list){

##1st

tmp<-GEP_reg_beta[,c("GENE", Program)]
colnames(tmp)<-c("GENE", "perturb_beta")
tmp<-data.frame(tmp, ensg=unname(symbol_to_ensg[as.character(tmp$GENE)]))
df<-merge(tmp, LOF, by="ensg")
df<-merge(df, shet, by="ensg")
df$post_mean<-scale(df$post_mean)
df$perturb_beta<-scale(df$perturb_beta)
df<-na.omit(df)

fit<- lm(post_mean~perturb_beta + shet, data=df)
P_withShet<-summary(fit)$coefficients[2,4]
beta_withShet<-summary(fit)$coefficients[2,1]
betaSE_withShet<-summary(fit)$coefficients[2,2]


hoge<-data.frame(Program=Program, cov="shet", 
    P_withShet=P_withShet, beta_withShet=beta_withShet, betaSE_withShet=betaSE_withShet)
summary_reg<-rbind(summary_reg, hoge)

##2nd
if(!is.element(Program, top_list[1])){
tmp<-GEP_reg_beta[,c("GENE", Program, top_list[1])]
colnames(tmp)<-c("GENE", "perturb_beta", "cov1")
tmp<-data.frame(tmp, ensg=unname(symbol_to_ensg[as.character(tmp$GENE)]))
df<-merge(tmp, LOF, by="ensg")
df<-merge(df, shet, by="ensg")
df$post_mean<-scale(df$post_mean)
df$perturb_beta<-scale(df$perturb_beta)
df<-na.omit(df)

fit<- lm(post_mean~perturb_beta + shet + cov1, data=df)
P_withShet<-summary(fit)$coefficients[2,4]
beta_withShet<-summary(fit)$coefficients[2,1]
betaSE_withShet<-summary(fit)$coefficients[2,2]

hoge<-data.frame(Program=Program, cov=paste(c("shet", top_list[1]), collapse=","), 
    P_withShet=P_withShet, beta_withShet=beta_withShet, betaSE_withShet=betaSE_withShet)
summary_reg<-rbind(summary_reg, hoge)
}


##3rd
if(!is.element(Program, top_list[1:2])){
tmp<-GEP_reg_beta[,c("GENE", Program, top_list[1:2])]
colnames(tmp)<-c("GENE", "perturb_beta", "cov1", "cov2")
tmp<-data.frame(tmp, ensg=unname(symbol_to_ensg[as.character(tmp$GENE)]))
df<-merge(tmp, LOF, by="ensg")
df<-merge(df, shet, by="ensg")
df$post_mean<-scale(df$post_mean)
df$perturb_beta<-scale(df$perturb_beta)
df<-na.omit(df)

fit<- lm(post_mean~perturb_beta + shet + cov1 + cov2, data=df)
P_withShet<-summary(fit)$coefficients[2,4]
beta_withShet<-summary(fit)$coefficients[2,1]
betaSE_withShet<-summary(fit)$coefficients[2,2]

hoge<-data.frame(Program=Program, cov=paste(c("shet", top_list[1:2]), collapse=","), 
    P_withShet=P_withShet, beta_withShet=beta_withShet, betaSE_withShet=betaSE_withShet)
summary_reg<-rbind(summary_reg, hoge)
}

##4th

if(!is.element(Program, top_list[1:3])){
tmp<-GEP_reg_beta[,c("GENE", Program, top_list[1:3])]
colnames(tmp)<-c("GENE", "perturb_beta", "cov1", "cov2", "cov3")
tmp<-data.frame(tmp, ensg=unname(symbol_to_ensg[as.character(tmp$GENE)]))
df<-merge(tmp, LOF, by="ensg")
df<-merge(df, shet, by="ensg")
df$post_mean<-scale(df$post_mean)
df$perturb_beta<-scale(df$perturb_beta)
df<-na.omit(df)

fit<- lm(post_mean~perturb_beta + shet + cov1 + cov2 + cov3, data=df)
P_withShet<-summary(fit)$coefficients[2,4]
beta_withShet<-summary(fit)$coefficients[2,1]
betaSE_withShet<-summary(fit)$coefficients[2,2]

hoge<-data.frame(Program=Program, cov=paste(c("shet", top_list[1:3]), collapse=","), 
    P_withShet=P_withShet, beta_withShet=beta_withShet, betaSE_withShet=betaSE_withShet)
summary_reg<-rbind(summary_reg, hoge)
}
}

write.table(summary_reg, paste0("multiCell_model/conditional/HbA1c_regulators_enrichment_top4.txt"), row.names=F, sep="\t", quote=F)



##make plot
##
df<-read.table("multiCell_model/Regulators_correlation.txt",header=T, stringsAsFactor=F)
##

library(reshape2)
mat<-dcast(df, A~B, value.var="R")
row.names(mat)<-mat[,1]
mat<-mat[,-1]
mat[is.na(mat)]<-1

library(pheatmap)
#pdf("regulators_correlation_heatmap.pdf")
p<-pheatmap(mat)
#p
#dev.off()

ORDER=row.names(mat)[p$tree_row$order]

##

data<-read.table(paste0("multiCell_model/conditional/HbA1c_regulators_enrichment_top4.txt"), header=T, stringsAsFactor=F)

top_list=unique(data$cov)

burden_1<-data[is.element(data$cov, top_list[1]),]

row.names(burden_1)<-burden_1[,1]
burden_1<-burden_1[colnames(mat),]

##
burden_2<-data[is.element(data$cov, top_list[2]),]

row.names(burden_2)<-burden_2[,1]
burden_2<-burden_2[colnames(mat),]
burden_2[is.na(burden_2)]<-1


##
burden_3<-data[is.element(data$cov, top_list[3]),]

row.names(burden_3)<-burden_3[,1]
#burden_1<-burden_1[rev(ORDER),]
burden_3<-burden_3[colnames(mat),]
burden_3[is.na(burden_3)]<-1

##
burden_4<-data[is.element(data$cov, top_list[4]),]

row.names(burden_4)<-burden_4[,1]
#burden_1<-burden_1[rev(ORDER),]
burden_4<-burden_4[colnames(mat),]
burden_4[is.na(burden_4)]<-1



##
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
col_fun = colorRamp2(c(-1, 0, 1 ), c("blue", "white", "red"))

tmp<-rbind(burden_1, burden_2, burden_3,burden_4)

tmp2<-data.frame(P=tmp$P_withShet, FDR=p.adjust(tmp$P_withShet, method="BH"))
tmp2<-tmp2[order(tmp2$P),]
THRESH=0.00106 ##set same value with the overall figure (FDR<0.01)
MAX=max(-log10(tmp$P_withShet))+5


burden_1<-data.frame(burden_1, COLOR=ifelse( (burden_1$P_withShet<THRESH & burden_1$beta_withShet >0), "red", ifelse ( (burden_1$P_withShet<THRESH & burden_1$beta_withShet <0), "blue", "grey")))
burden_2<-data.frame(burden_2, COLOR=ifelse( (burden_2$P_withShet<THRESH & burden_2$beta_withShet >0), "red", ifelse ( (burden_2$P_withShet<THRESH & burden_2$beta_withShet <0), "blue", "grey")))
burden_3<-data.frame(burden_3, COLOR=ifelse( (burden_3$P_withShet<THRESH & burden_3$beta_withShet >0), "red", ifelse ( (burden_3$P_withShet<THRESH & burden_3$beta_withShet <0), "blue", "grey")))
burden_4<-data.frame(burden_4, COLOR=ifelse( (burden_4$P_withShet<THRESH & burden_4$beta_withShet >0), "red", ifelse ( (burden_4$P_withShet<THRESH & burden_4$beta_withShet <0), "blue", "grey")))


cell_type=ifelse(grepl("K562GW", colnames(mat)), "K562GW", ifelse(grepl("K562_essential", colnames(mat)), "K562essential", ifelse(grepl("hepg2", colnames(mat)), "HepG2", ifelse(grepl("jurkat", colnames(mat)), "Jurkat", 
    ifelse(grepl("GSM6858449", colnames(mat)),"THP1", "RPE1")))))

column_ha3 = HeatmapAnnotation(
  "All correlations" = anno_barplot(-log10(burden_1$P_withShet),  gp = gpar(fill = burden_1$COLOR), ylim=c(0,MAX), axis_param=list(gp=gpar(fontsize = 20))), 
  "top1 regressed" = anno_barplot(-log10(burden_2$P_withShet),  gp = gpar(fill = burden_2$COLOR), ylim=c(0,MAX), axis_param=list(gp=gpar(fontsize = 20))), 
  "top2 regressed" = anno_barplot(-log10(burden_3$P_withShet),  gp = gpar(fill = burden_3$COLOR), ylim=c(0,MAX), axis_param=list(gp=gpar(fontsize = 20))), 
  "top3 regressed" = anno_barplot(-log10(burden_4$P_withShet),  gp = gpar(fill = burden_4$COLOR), ylim=c(0,MAX), axis_param=list(gp=gpar(fontsize = 20))), 
  "K562GW" = ifelse(cell_type=="K562GW", 1, 0), "K562essential" = ifelse(cell_type=="K562essential", 1, 0),
   "HepG2" = ifelse(cell_type=="HepG2", 1, 0),
    "Jurkat" = ifelse(cell_type=="Jurkat", 1, 0),
     "RPE1" = ifelse(cell_type=="RPE1", 1, 0),
  col = list("K562GW" = c("1" = "black", "0"="white"), "K562essential" = c("1" = "black", "0"="white"), "HepG2" = c("1" = "black", "0"="white"), "Jurkat" = c("1" = "black", "0"="white"), "RPE1" = c("1" = "black", "0"="white")),
height =  unit(200, "mm"), show_legend = FALSE, annotation_name_gp= gpar(fontsize = 20))

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
         top_annotation = column_ha3)
##

##
df2<-read.table("multiCell_model/Program_genes_Jaccard.txt",header=T, stringsAsFactor=F)

library(reshape2)
mat2<-dcast(df2, A~B, value.var="JI")
row.names(mat2)<-mat2[,1]
mat2<-mat2[,-1]
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
        }
        },
        heatmap_legend_param = list(
        title = "Program\nJaccard index", 
        labels_gp=gpar(fontsize = 15), title_gp=gpar(fontsize = 15)
        )
         #rect_gp = gpar(col = "black", lwd = 1),　
         #top_annotation = column_ha3
         )



pdf("multiCell_model/conditional/Combined_HbA1c_conditionalPlot_forSupplement.pdf", height=13, width=15)
set.seed(50)
draw(ht1 + ht2, ht_gap = unit(-297, "mm") ,padding = unit(c(0,5,0, 35), "mm"))
dev.off()


