library(limma)
library(ggplot2)
library(ggpubr)

gene="UCHL5"                
expFile="exp.txt"         
clikFile="clinical.txt"    
setwd("File path")      

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

qx=as.numeric(quantile(data, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ((qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
	data[data<0]=0
	data=log2(data+1)}
#data=normalizeBetweenArrays(data)

exp=t(data[gene,,drop=F])

clinical=read.table(clikFile, header=T, sep="\t", check.names=F, row.names=1)
cliName=colnames(clinical)[1]
	
sameSample=intersect(row.names(exp), row.names(clinical))
exp=exp[sameSample,,drop=F]
clinical=clinical[sameSample,,drop=F]
data=cbind(as.data.frame(exp), as.data.frame(clinical))

group=levels(factor(data[,cliName]))
data[,cliName]=factor(data[,cliName], levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

boxplot=ggboxplot(data, x=cliName, y=gene, fill=cliName,
			      xlab="GSE67501",
			      ylab=paste0(gene, " expression"),
			      legend.title="",
			      palette = c("blue", "red") )+ 
	stat_compare_means(comparisons = my_comparisons)	

pdf(file=paste0(gene, ".boxplot.pdf"), width=5, height=4.5)
print(boxplot)
dev.off()