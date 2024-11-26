library(fmsb) 
expFile="UCHL5geneExp.txt"      
tmbFile="TMB.txt"         
col="red"                  
setwd("File path")      

exp=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
exp=exp[(exp[,"Type"]=="Tumor"),]

TMB=read.table(tmbFile, header=T, sep="\t", check.names=F, row.names=1)
TMB=as.matrix(TMB)
row.names(TMB)=gsub(".$", "", row.names(TMB))

sameSample=intersect(row.names(TMB), row.names(exp))
TMB=TMB[sameSample,]
exp=exp[sameSample,]

outTab=data.frame()
fmsbTab=data.frame()
for(i in levels(factor(exp[,"CancerType"]))){
     exp1=exp[(exp[,"CancerType"]==i),]
     TMB1=TMB[(TMB[,"CancerType"]==i),]
	 x=as.numeric(TMB1[,1])
	 y=as.numeric(exp1[,1])
	 corT=cor.test(x,y,method="spearman")
	 cor=corT$estimate
	 pValue=corT$p.value
	 sig=ifelse(pValue<0.001,"***",ifelse(pValue<0.01,"**",ifelse(pValue<0.05,"*"," ")))
	 outTab=rbind(outTab, cbind(CancerType=i,cor=cor,pValue=pValue))
	 fmsbTab=rbind(fmsbTab, cbind(CancerType=paste0(i,sig),cor=cor))
}
write.table(outTab,file="cor.result.txt",sep="\t",row.names=F,quote=F)
write.table(t(fmsbTab),file="fmsbInput.txt",sep="\t",col.names=F,quote=F)

data=read.table("fmsbInput.txt", header=T, sep="\t", check.names=F, row.names=1)
maxValue=ceiling(max(abs(data))*10)/10
data=rbind(rep(maxValue,ncol(data)),rep(-maxValue,ncol(data)),data)

pdf(file="radar.pdf", width=7, height=7)
radarchart(data, axistype=1,
	title="Tumor Mutation Burden",   
    pcol=col,                   
    plwd=2 ,                     
    plty=1,                      
    cglcol="grey",               
    cglty=1,                     
    caxislabels=seq(-maxValue,maxValue,maxValue/2),   
    cglwd=1.2,                  
    axislabcol="blue",           
    vlcex=0.8                    
)
dev.off()