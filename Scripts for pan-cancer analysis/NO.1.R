library(limma)    
gene="UCHL5"        
setwd("File path")    

files=dir()
files=grep("^symbol.", files, value=T)

geneList=list()
for(i in files){
	CancerType=gsub("symbol\\.|\\.txt", "", i)
	rt=read.table(i, header=T, sep="\t", check.names=F)
    geneList[[CancerType]]=as.vector(rt[,1])
}
interGenes=Reduce(intersect, geneList)

outTab=data.frame()
allTab=data.frame()
j=0
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
	Type=ifelse(group==0, "Tumor", "Normal")
	
	geneExp=t(data[gene,,drop=F])
	outTab=rbind(outTab, cbind(geneExp,Type,CancerType))
	
	j=j+1
    if(j==1){
    	allTab=data[interGenes,group==0]
    }else{
    	allTab=cbind(allTab, data[interGenes,group==0])
    }
}

out=cbind(ID=row.names(outTab), outTab)
write.table(out, file="geneExp.txt", sep="\t", quote=F, row.names=F)

allTab=allTab[rowMeans(allTab)>0.1,]
allTab=cbind(ID=row.names(allTab), allTab)
write.table(allTab, file="merge.txt", sep="\t", quote=F, row.names=F)