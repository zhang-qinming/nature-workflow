args   <- commandArgs(trailingOnly = T)   # enable reading factors from linux command line
options("scipen"=10) ##disable 100000 to 1e6
TOP_N<-as.character( args[1] ) ##100,200,300

TRAIT=if (length(args) >= 2) as.character(args[2]) else "MCH"
example_file <- if (length(args) >= 3) as.character(args[3]) else paste0(TRAIT, "_program5_regulator3_LOF0.1.txt")

summary<-data.frame()
file_list=dir(paste0("permutation_test/P", TOP_N))
file_list=file_list[grep(paste0(TRAIT), file_list)]
file_list=file_list[-grep("Selected", file_list)]

for(FILE in file_list){
cat(paste0("START ", FILE, "\n"))
data<-read.table(paste0("permutation_test/P", TOP_N, "/", FILE), header=T, stringsAsFactor=F)

bin1_concoN=data[1,3]
bin1_fisherP=fisher.test(matrix(c(data[1,3], data[1,5]-data[1,3], data[1,6], data[1,8]-data[1,6]), nrow=2, byrow=T))$p.value
random_bin1<-c()
seed_list=unique(data$SEED)
seed_list=seed_list[!is.element(seed_list, "TRUE")]

if(length(seed_list) > 19995){
for(S in seed_list){
tmp<-data[is.element(data$SEED, S),]
tmp<-tmp[is.element(tmp$BIN, 1),,drop=F]
tmp<-na.omit(tmp)
if(ncol(tmp)==8 && !is.na(tmp[1,8])){
if(tmp[1,8]-tmp[1,6] >0){
tmp_fisherP=fisher.test(matrix(c(tmp[1,3], tmp[1,5]-tmp[1,3], tmp[1,6], tmp[1,8]-tmp[1,6]), nrow=2, byrow=T))$p.value
random_bin1<-c(random_bin1, tmp_fisherP)
}}}

if(length(random_bin1) > 19995){
P=rank(c(bin1_fisherP, random_bin1))[1]/length(random_bin1)

hoge<-data.frame(FILE=FILE, TRAIT=TRAIT, BIN="bin1", P=P, N=bin1_concoN, PermN=length(random_bin1))
summary<-rbind(summary, hoge)
}
}
}
##
##
write.table(summary, paste0("P", TOP_N, "_", TRAIT, "_permP.sum"), row.names=F, sep="\t", quote=F)


##make figure

data<-read.table(paste0("permutation_test/P", TOP_N, "/", example_file), header=T, stringsAsFactor=F)
data<-data[is.element(data$BIN, 1),]

summary<-data.frame()
bin1_concoN=data[1,3]
bin1_fisherP=fisher.test(matrix(c(data[1,3], data[1,5]-data[1,3], data[1,6], data[1,8]-data[1,6]), nrow=2, byrow=T))$p.value
hoge<-data.frame(SEED="TRUE", fisher_P=bin1_fisherP)
summary<-rbind(summary, hoge)
seed_list=unique(data$SEED)
seed_list=seed_list[!is.element(seed_list, "TRUE")]
for(S in seed_list){
tmp<-data[is.element(data$SEED, S),]
tmp<-tmp[is.element(tmp$BIN, 1),,drop=F]
tmp_fisherP=fisher.test(matrix(c(tmp[1,3], tmp[1,5]-tmp[1,3], tmp[1,6], tmp[1,8]-tmp[1,6]), nrow=2, byrow=T))$p.value
hoge<-data.frame(SEED=S, fisher_P=tmp_fisherP)
summary<-rbind(summary, hoge)
}

summary<-data.frame(summary, logP=-log10(summary$fisher_P))
library(ggplot2)
library(ggarchery)
g<-ggplot(summary, aes(x=logP))
g<-g +  theme_classic(base_size = 30, base_family = "Helvetica")
g<-g + geom_histogram(fill="black", position="dodge", binwidth=0.05)
g<-g + geom_arrowsegment(aes(x = summary$logP[1], xend = summary$logP[1], y = 100, yend = 0), arrows = arrow(type = "closed")) 
g<-g + ylab("count") + xlab("log10(Fisher P)")
ggsave(paste0("P", TOP_N, "_", sub(".txt$", "", example_file), "_hist.pdf"))

##

data_f<-read.table(paste0("permutation_test/P", TOP_N, "/", example_file), header=T, stringsAsFactor=F)
data_f<-data.frame(data_f, COLOR=ifelse(data_f$SEED=="TRUE", 1, 0))
data_ff<-data_f[is.element(data_f$COLOR, 1),]

library(ggplot2)
g<-ggplot(data_f, aes(x=CONCO_N, y=DISCO_N, color=factor(COLOR)))
g<-g +  theme_classic(base_size = 30, base_family = "Helvetica")
g<-g + geom_point(size=1, alpha=0.5)
g<-g + geom_point(data=data_ff, size=5, color="red")
g<-g + scale_color_manual(values=c("grey", "red"))
g<-g + xlab("N, concordant genes") + ylab("N, discordant genes")
g<-g + stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white")
g<-g + theme(legend.position="none")
ggsave(plot=g, paste0(sub(".txt$", "", example_file), "_concordant_discordant_gene_num_scatter.pdf"), height=8, width=8)
