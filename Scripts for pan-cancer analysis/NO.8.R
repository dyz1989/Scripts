library(limma)      
pFilter=0.05        
setwd("File path")       

files=dir()
files=grep("^symbol.", files, value=T)

outTab=data.frame()
for(i in files){
	CancerType=gsub("symbol\\.|\\.txt", "", i)
	rt=read.table(i, header=T, sep="\t", check.names=F)
	
	rt=as.matrix(rt)
	rownames(rt)=rt[,1]
	exp=rt[,2:ncol(rt)]
	dimnames=list(rownames(exp), colnames(exp))
	data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
	data=t(avereps(data))
	row.names(data)=gsub(".$", "", row.names(data))
	data=t(avereps(data))

	group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
	group=sapply(strsplit(group,""), "[", 1)
	group=gsub("2", "1", group)
	data=data[,group==0]
	data=data[rowMeans(data)>0,]

	v <-voom(data, plot=F, save.plot=F)
	out=v$E
	out=rbind(ID=colnames(out),out)
	write.table(out, file="uniq.symbol.txt", sep="\t", quote=F, col.names=F)      
	
	source("panImmune15.CIBERSORT.R")
	results=CIBERSORT("ref.txt", "uniq.symbol.txt", perm=1000, QN=TRUE)

	immune=read.table("CIBERSORT-Results.txt", header=T, sep="\t", check.names=F, row.names=1)
	immune=immune[immune[,"P-value"]<pFilter,]
	immune=as.matrix(immune[,1:(ncol(immune)-3)])
	outTab=rbind(outTab, cbind(immune,CancerType))
	file.remove("CIBERSORT-Results.txt")
	file.remove("uniq.symbol.txt")
}
out=cbind(ID=row.names(outTab), outTab)
write.table(out, file="CIBERSORT.result.txt", sep="\t", quote=F, row.names=F)