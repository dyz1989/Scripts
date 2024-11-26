library(limma)            
expFile="UCHL5geneExp.txt"     
setwd("File path")    

files=dir()
files=grep(".survival.tsv$", files, value=T)

surTime=data.frame()
for(surFile in files){
    rt=read.table(surFile, header=T, sep="\t", check.names=F, row.names=1)
    rt=rt[,c(3,1)]
    surTime=rbind(surTime, rt)
}
colnames(surTime)=c("futime","fustat")
surTime=as.matrix(surTime)
row.names(surTime)=gsub(".$","",row.names(surTime))

exp=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
exp=exp[(exp[,"Type"]=="Tumor"),]

sameSample=intersect(row.names(surTime), row.names(exp))
surTime=surTime[sameSample,,drop=F]
exp=exp[sameSample,,drop=F]
surData=cbind(surTime, exp)
surData=cbind(id=row.names(surData), surData)
write.table(surData, file="expTime.txt", sep="\t", quote=F, row.names=F)



