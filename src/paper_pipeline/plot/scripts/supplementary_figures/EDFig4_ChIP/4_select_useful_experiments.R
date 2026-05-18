args   <- commandArgs(trailingOnly = T)  

##For each TF, choose the best metric for predicting the TF effect.

res<-read.table("data/TF_ChIP/chip_score_KOeffects_MWUtest.txt", header=T, stringsAsFacto=F)

library(dplyr)
res_f<-res %>% filter(SCORE=="chip_score")   %>% filter( G==1000 | G==5000 | G==10000 | G==50000 | G==100000 | G==500000 | G==1000000 ) 

P_list<-sort(na.omit(c(res_f$MWU_up,res_f$MWU_dn )))
thresh<-P_list[p.adjust(P_list, method="BH")>0.05][1]

final<-data.frame()
TF_list=unique(res_f$TF)

for(TF in TF_list){
tmp<-res_f[is.element(res_f$TF, TF),]
if(length(na.omit(c(tmp$MWU_up,tmp$MWU_dn )))>0){
minP<-min(na.omit(c(tmp$MWU_up,tmp$MWU_dn )))

if(minP<thresh){
tmp2<- tmp %>% filter(tmp$MWU_up== minP | tmp$MWU_dn== minP ) 
TYPE=ifelse(tmp2$MWU_up[1] < thresh & tmp2$MWU_dn[1] < thresh & !is.na(tmp2$MWU_dn[1]) & !is.na(tmp2$MWU_up[1]) , "both", 
	ifelse(tmp2$MWU_up[1] < thresh & !is.na(tmp2$MWU_up[1]), "repressor", "activator"))
tmp2<-data.frame(tmp2, TYPE=TYPE)
final<-rbind(final, tmp2)
}
}}

final<-final[order(final$TYPE),]
write.table(final, "data/TF_ChIP/selected_scores_for_prediction.txt", row.names=F, sep="\t", quote=F)
