library(ggplot2)
library(ggpubr)
library(ggExtra)

corFilter=0.5         
pFilter=0.001        
expFile="UCHL5geneExp.txt"            
immFile="CIBERSORT.result.txt"    
setwd("File path")      

exp=read.table(expFile, header=T,sep="\t", check.names=F, row.names=1)
exp=exp[(exp[,"Type"]=="Tumor"),]
gene=colnames(exp)[1]

immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)

sameSample=intersect(row.names(immune), row.names(exp))
immune=immune[sameSample,]
exp=exp[sameSample,]

outTab=data.frame()
for(i in levels(factor(exp[,"CancerType"]))){
    exp1=exp[(exp[,"CancerType"]==i),]
    immune1=immune[(immune[,"CancerType"]==i),]
    y=as.numeric(exp1[,1])
    outVector=data.frame(i, gene)
	for(j in colnames(immune1)[1:22]){
		x=as.numeric(immune1[,j])
		if(sd(x)>0.01){
			df1=as.data.frame(cbind(x,y))
			corT=cor.test(x,y,method="spearman")
			cor=corT$estimate
			pValue=corT$p.value
			outVector=cbind(outVector,pValue)
			
			if((abs(cor)>corFilter) & (pValue<pFilter)){
				p1=ggplot(df1, aes(x, y)) + 
					xlab(j)+ylab(gene)+
					ggtitle(paste0("\nCancer: ",i))+theme(title=element_text(size=10))+
				    geom_point()+ geom_smooth(method="lm", formula=y ~ x) + theme_bw()+
				    stat_cor(method = 'spearman', aes(x =x, y =y))
			    p2=ggMarginal(p1, type = "density", xparams = list(fill = "orange"),yparams = list(fill = "blue"))
				pdf(file=paste0("estimateCor.",i,"_",j,".pdf"), width=5, height=5.1)
				print(p2)
				dev.off()
			}
		}
		else{
			outVector=cbind(outVector,pValue=1)
		}
	}
	outTab=rbind(outTab, outVector)
}

colNames=c("CancerType", "Gene", colnames(immune)[1:22])
colnames(outTab)=colNames
write.table(outTab, file="immuneCor.result.txt", sep="\t", row.names=F, quote=F)