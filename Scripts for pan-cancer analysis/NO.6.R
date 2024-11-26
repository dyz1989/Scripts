library(limma)
library(estimate)
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
	dimnames=list(rownames(exp),colnames(exp))
	data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
	data=t(avereps(data))
	row.names(data)=gsub(".$", "", row.names(data))
	data=t(avereps(data))
	group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
	group=sapply(strsplit(group,""), "[", 1)
	group=gsub("2", "1", group)
	data=data[,group==0]
	out=data[rowMeans(data)>0,]
	out=rbind(ID=colnames(out), out)
	write.table(out, file="uniq.symbol.txt", sep="\t", quote=F, col.names=F)
	
	filterCommonGenes(input.f="uniq.symbol.txt", output.f="commonGenes.gct", id="GeneSymbol")
	estimateScore(input.ds="commonGenes.gct", output.ds="estimateScore.gct")
	
	scores=read.table("estimateScore.gct", header=T, skip=2)
	rownames(scores)=scores[,1]
	scores=t(scores[,3:ncol(scores)])
	rownames(scores)=gsub("\\.", "\\-", rownames(scores))
	outTab=rbind(outTab,cbind(scores, CancerType))
	file.remove("uniq.symbol.txt")
	file.remove("commonGenes.gct")
	file.remove("estimateScore.gct")
}
out=cbind(ID=row.names(outTab), outTab)
write.table(out, file="estimateScores.txt", sep="\t", quote=F, row.names=F)