args   <- commandArgs(trailingOnly = T)   # enable reading factors from linux command line
trait <- as.character( args[1] )

TOP<-read.table(paste0("clump/", trait, ".top.100kbpmerged"), header=T, stringsAsFactor=F)
TOP<-TOP[order(TOP$POS),]
TOP<-TOP[order(TOP$CHR),]

TOP<-data.frame(CHR=paste0("chr", TOP$CHR), POS1=TOP$POS, POS2=TOP$POS + 1, P=TOP$P )

write.table(TOP, paste0("tmp_", trait, "_top.bed"), row.names=F, sep="\t", quote=F, col.names=F)