args   <- commandArgs(trailingOnly = T)  
K<-as.numeric( args[1] ) ##60

library(data.table)
df<-fread("data/Perturbseq/metadata/gwps_metadata.csv", header=T)


##data loading
GEP<-read.table(paste0("data/Perturbseq/cNMF/K562GW/cNMF_all.gene_spectra_score.k_", K, ".dt_0_5.txt"), header=T, stringsAsFactor=F)
GEP<-t(GEP)
colnames(GEP)<-paste0("P", c(1:K))
corresp<-read.table("data/gencode_v41_gname_gid_ALL_sorted_onlyID", header=F, stringsAsFactor=F)
corresp<-corresp[!duplicated(corresp[,1]),]
row.names(corresp)<-corresp[,1]
GEP<-data.frame(GENE=corresp[row.names(GEP), 2], GEP)

GEP_reg_beta<-data.frame()
GEP_reg_P<-data.frame()

for(i in 1:K){
tmp<-read.table(paste0("data/Perturbseq/cNMF_regulation/K562GW/K",  K, "_program", i, "_perturb_effects.txt"), header=T)

tmp1<-tmp[,c(1,2)]
tmp2<-tmp[,c(1,3)]

colnames(tmp1)<-c("GENE", paste0("P", i))
colnames(tmp2)<-c("GENE", paste0("P", i))

if(i==1){
GEP_reg_beta<-tmp1
GEP_reg_P<-tmp2
} else {
GEP_reg_beta<-merge(GEP_reg_beta, tmp1, by="GENE")
GEP_reg_P<-merge(GEP_reg_P, tmp2, by="GENE")
}
}

##
corresp<-read.table("data/gencode_v41_gname_gid_ALL_sorted_onlyID", header=F, stringsAsFactor=F)
corresp<-corresp[!duplicated(corresp[,1]),]
row.names(corresp)<-corresp[,1]

summary<-data.frame()
T1="P25"
target<-unique(c(GEP_reg_P$GENE[p.adjust(GEP_reg_P[,T1], method="BH")<0.05]))
library(dplyr)
Ra_genes<-GEP_reg_beta[is.element(GEP_reg_beta$GENE, target),] %>% filter(P25 < 0) %>% select(GENE)
Rb_genes<-GEP_reg_beta[is.element(GEP_reg_beta$GENE, target),]  %>% filter(P25 > 0) %>% select(GENE) 

P4_program<-GEP[order(GEP[,"P4"], decreasing=T), "GENE"][1:200]
P25_program<-GEP[order(GEP[,"P25"], decreasing=T), "GENE"][1:200]

df_Ra<-df %>% filter(is.element(gene, Ra_genes[,1]))
df_Rb<-df %>% filter(is.element(gene, Rb_genes[,1]))

df_P4<-df %>% filter(is.element(gene, P4_program))
df_P25<-df %>% filter(is.element(gene, P25_program))

df_ctrl<-df %>% filter(is.element(gene, "non-targeting"))

Ra_res<-data.frame(Program="G2Mphase", GROUP="up_regulator", table(df_Ra$phase)/sum(table(df_Ra$phase)))
Rb_res<-data.frame(Program="G2Mphase", GROUP="dn_regulator", table(df_Rb$phase)/sum(table(df_Rb$phase)))

summary<-rbind(summary, Ra_res)
summary<-rbind(summary, Rb_res)

##

T1="P4"
target<-unique(c(GEP_reg_P$GENE[p.adjust(GEP_reg_P[,T1], method="BH")<0.05]))
library(dplyr)
Ra_genes<-GEP_reg_beta[is.element(GEP_reg_beta$GENE, target),] %>% filter(P4 < 0) %>% select(GENE)
Rb_genes<-GEP_reg_beta[is.element(GEP_reg_beta$GENE, target),]  %>% filter(P4 > 0) %>% select(GENE) 

P4_program<-GEP[order(GEP[,"P4"], decreasing=T), "GENE"][1:200]
P25_program<-GEP[order(GEP[,"P25"], decreasing=T), "GENE"][1:200]

df_Ra<-df %>% filter(is.element(gene, Ra_genes[,1]))
df_Rb<-df %>% filter(is.element(gene, Rb_genes[,1]))

df_P4<-df %>% filter(is.element(gene, P4_program))
df_P25<-df %>% filter(is.element(gene, P25_program))

df_ctrl<-df %>% filter(is.element(gene, "non-targeting"))

Ra_res<-data.frame(Program="Sphase", GROUP="up_regulator", table(df_Ra$phase)/sum(table(df_Ra$phase)))
Rb_res<-data.frame(Program="Sphase", GROUP="dn_regulator", table(df_Rb$phase)/sum(table(df_Rb$phase)))

summary<-rbind(summary, Ra_res)
summary<-rbind(summary, Rb_res)

P4_res<-data.frame(Program="Sphase", GROUP="Programgenes", table(df_P4$phase)/sum(table(df_P4$phase)))
P25_res<-data.frame(Program="G2Mphase", GROUP="Programgenes", table(df_P25$phase)/sum(table(df_P25$phase)))
ctrl_res<-data.frame(Program="ctrl", GROUP="non-targeting",table(df_ctrl$phase)/sum(table(df_ctrl$phase)))

summary<-rbind(summary, P4_res)
summary<-rbind(summary, P25_res)
summary<-rbind(summary, ctrl_res)



##calculate jackknife variance


summary2<-data.frame()
T1="P25"
target<-unique(c(GEP_reg_P$GENE[p.adjust(GEP_reg_P[,T1], method="BH")<0.05]))
library(dplyr)
Ra_genes<-GEP_reg_beta[is.element(GEP_reg_beta$GENE, target),] %>% filter(P25 < 0) %>% select(GENE)
Rb_genes<-GEP_reg_beta[is.element(GEP_reg_beta$GENE, target),]  %>% filter(P25 > 0) %>% select(GENE) 

P4_program<-GEP[order(GEP[,"P4"], decreasing=T), "GENE"][1:200]
P25_program<-GEP[order(GEP[,"P25"], decreasing=T), "GENE"][1:200]

df_ctrl<-df %>% filter(is.element(gene, "non-targeting"))
ctrl_res<-data.frame(Program="ctrl", GROUP="non-targeting",GENE="non-targeting", table(df_ctrl$phase)/sum(table(df_ctrl$phase)))

summary2<-rbind(summary2, ctrl_res)

for(GENE in Ra_genes[,1]){
df_Ra<-df %>% filter(is.element(gene, GENE))
if(nrow(df_Ra)>30){
Ra_res<-data.frame(Program="G2Mphase", GROUP="up_regulator", GENE=GENE, table(df_Ra$phase)/sum(table(df_Ra$phase)))
summary2<-rbind(summary2, Ra_res)
}}

for(GENE in Rb_genes[,1]){
df_Ra<-df %>% filter(is.element(gene, GENE))
if(nrow(df_Ra)>30){
Ra_res<-data.frame(Program="G2Mphase", GROUP="dn_regulator", GENE=GENE, table(df_Ra$phase)/sum(table(df_Ra$phase)))
summary2<-rbind(summary2, Ra_res)
}}

for(GENE in P4_program){
df_Ra<-df %>% filter(is.element(gene, GENE))
if(nrow(df_Ra)>30){
Ra_res<-data.frame(Program="Sphase", GROUP="Programgenes", GENE=GENE, table(df_Ra$phase)/sum(table(df_Ra$phase)))
summary2<-rbind(summary2, Ra_res)
}}

for(GENE in P25_program){
df_Ra<-df %>% filter(is.element(gene, GENE))
if(nrow(df_Ra)>30){
Ra_res<-data.frame(Program="G2Mphase", GROUP="Programgenes", GENE=GENE, table(df_Ra$phase)/sum(table(df_Ra$phase)))
summary2<-rbind(summary2, Ra_res)
}}


##

T1="P4"
target<-unique(c(GEP_reg_P$GENE[p.adjust(GEP_reg_P[,T1], method="BH")<0.05]))
library(dplyr)
Ra_genes<-GEP_reg_beta[is.element(GEP_reg_beta$GENE, target),] %>% filter(P4 < 0) %>% select(GENE)
Rb_genes<-GEP_reg_beta[is.element(GEP_reg_beta$GENE, target),]  %>% filter(P4 > 0) %>% select(GENE) 


for(GENE in Ra_genes[,1]){
df_Ra<-df %>% filter(is.element(gene, GENE))
if(nrow(df_Ra)>30){
Ra_res<-data.frame(Program="Sphase", GROUP="up_regulator", GENE=GENE, table(df_Ra$phase)/sum(table(df_Ra$phase)))
summary2<-rbind(summary2, Ra_res)
}}

for(GENE in Rb_genes[,1]){
df_Ra<-df %>% filter(is.element(gene, GENE))
if(nrow(df_Ra)>30){
Ra_res<-data.frame(Program="Sphase", GROUP="dn_regulator", GENE=GENE, table(df_Ra$phase)/sum(table(df_Ra$phase)))
summary2<-rbind(summary2, Ra_res)
}}



##JackKnife_variance

S_var=0
G2M_var=0
G1_var=0


summary3<-data.frame()
for(P in c("Sphase", "G2Mphase")){
for(G in c("up_regulator", "dn_regulator", "Programgenes")){
df<-summary2  %>% filter(is.element(Program, P)) %>% filter(is.element(GROUP, G))
obs_S=apply(df %>% filter(is.element(Var1, "S"))  %>% select(Freq), 2, mean)
obs_G2M=apply(df %>% filter(is.element(Var1, "G2M"))  %>% select(Freq), 2, mean)
obs_G1=apply(df %>% filter(is.element(Var1, "G1"))  %>% select(Freq), 2, mean)

for(i in 1:nrow(df)){
 tmp<-df[-i,]
 tmp_S<- apply(tmp %>% filter(is.element(Var1, "S"))  %>% select(Freq), 2, mean)
 tmp_G2M<-apply(tmp %>% filter(is.element(Var1, "G2M"))  %>% select(Freq), 2, mean)
 tmp_G1<-apply(tmp %>% filter(is.element(Var1, "G1"))  %>% select(Freq), 2, mean)

 S_var=S_var + (tmp_S - obs_S)^2
 G2M_var=G2M_var + (tmp_G2M - obs_G2M)^2
 G1_var=G1_var + (tmp_G1 - obs_G1)^2

}
S_se=sqrt(S_var)
G2M_se=sqrt(G2M_var)
G1_se=sqrt(G1_var)

hoge<-data.frame(Program=P, Group=G, Sfreq=obs_S, G2Mfreq=obs_G2M, G1freq=obs_G1,  Sfreq_se=S_se, G2Mfreq_se=G2M_se, G1freq_se=G1_se )
summary3<-rbind(summary3, hoge)
##
}}



##make bar plot
mat<-summary
mat$Var1<-factor(as.character(mat$Var1), levels=c("G2M", "G1", "S"))
mat$Program<-factor(as.character(mat$Program), levels=c("ctrl", "Sphase", "G2Mphase"))

tmp<-mat[grep("ctrl", mat$Program),]
y_S=tmp[is.element(tmp$Var1, "S"),"Freq"]
y_M=tmp[is.element(tmp$Var1, "S"),"Freq"] + tmp[is.element(tmp$Var1, "G1"),"Freq"]

mat2<-summary3
mat2<-data.frame(mat2, comp=paste0(mat2$Program, "_", mat2$Group), Freq=1)
for(i in 1:nrow(mat2)){
P=mat2[i,1]
G=mat2[i,2]
S= as.numeric(mat  %>% filter(is.element(Program, P))  %>% filter(is.element(GROUP, G))  %>% filter(is.element(Var1, "S"))   %>% select(Freq) )
G1= as.numeric(mat  %>% filter(is.element(Program, P))  %>% filter(is.element(GROUP, G))  %>% filter(is.element(Var1, "G1"))   %>% select(Freq) )
G2M= as.numeric(mat  %>% filter(is.element(Program, P))  %>% filter(is.element(GROUP, G))  %>% filter(is.element(Var1, "G2M"))   %>% select(Freq) )
mat2[i,3]<-S
mat2[i,4]<-G2M
mat2[i,5]<-G1
}

library(ggplot2)
library(ggsci)
mat<-data.frame(mat, comp=paste0(mat$Program, "_", mat$GROUP))
mat$comp<-factor(as.character(mat$comp), levels=c("Sphase_Programgenes", "G2Mphase_Programgenes", "Sphase_up_regulator", "Sphase_dn_regulator",  "G2Mphase_up_regulator", "G2Mphase_dn_regulator",  "ctrl_non-targeting") )

g<-ggplot(mat, aes( y=comp, x=Freq))
g<-g +  theme_classic(base_size = 40, base_family = "Helvetica")
g<-g + geom_bar(aes(fill=Var1), position="stack", stat="identity")
g<-g + geom_vline(xintercept=y_S, linetype="dashed")
g<-g + geom_vline(xintercept=y_M, linetype="dashed")
g<-g + scale_fill_manual(values=c("S"="#FEA601", "G1"="#4783B5", "G2M"="#34a853"))
g<-g + geom_errorbarh(data=mat2, aes(xmin= Sfreq-  Sfreq_se, 
                    xmax=  Sfreq+ Sfreq_se), height=0, linewidth=1,
                position="identity")
g<-g + geom_errorbarh(data=mat2, aes(xmin= 1 - G2Mfreq- G2Mfreq_se, 
                    xmax=  1 - G2Mfreq+ G2Mfreq_se), height=0,linewidth=1,
                position="identity")
g<-g + ylab("") + xlab("Frequency")
ggsave(plot=g, "Fig5I.pdf", width=15, height=8)
