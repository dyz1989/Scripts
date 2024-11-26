library(survival)
library(survminer)
library(forestplot)

inputFile="expTime.txt"       
setwd("File path")       
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)    
rt$futime=rt$futime/365
gene=colnames(rt)[3]

outTab=data.frame()
for(i in levels(factor(rt[,"CancerType"]))){
	rt1=rt[(rt[,"CancerType"]==i),]
	
	cox=coxph(Surv(futime, fustat) ~ rt1[,gene], data = rt1)
	coxSummary = summary(cox)
    coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	outTab=rbind(outTab,
	             cbind(cancer=i,
	                   HR=coxSummary$conf.int[,"exp(coef)"],
	                   HR.95L=coxSummary$conf.int[,"lower .95"],
	                   HR.95H=coxSummary$conf.int[,"upper .95"],
			           pvalue=coxP) )
	
	group=ifelse(rt1[,gene]>median(rt1[,gene]), "high", "low")
	diff=survdiff(Surv(futime, fustat) ~ group, data=rt1)
	pValue=1-pchisq(diff$chisq,df=1)
	if(pValue<0.05){
		if(pValue<0.001){
			pValue="p<0.001"
		}else{
			pValue=paste0("p=",sprintf("%.03f",pValue))
		}
		
		fit=survfit(Surv(futime, fustat) ~ group, data = rt1)
		surPlot=ggsurvplot(fit, 
				    data=rt1,
				    title=paste0("Cancer: ",i),
				    pval=pValue,
				    pval.size=6,
				    conf.int=F,
				    legend.title=paste0(gene," levels"),
				    legend.labs=c("high","low"),
				    font.legend=12,
				    fontsize=4,
				    xlab="Time(years)",
				    ylab="Overall survival",
				    break.time.by = 2,
				    palette=c("#6495ED","#FF8C00"),
				    risk.table=TRUE,
				    risk.table.title="",
				    risk.table.height=.25)
		pdf(file=paste0("survival.",i,".pdf"), width=6, height=5, onefile=FALSE)
		print(surPlot)
		dev.off()
	}
}
write.table(outTab, file="cox.result.txt", sep="\t", row.names=F, quote=F)

inputFile="cox.result.txt"        
outFile="forest.pdf"         
setwd("File path")    

library(fdrtool)

rt=read.table(inputFile,header=T,sep="\t",row.names=1,check.names=F)
rt$qvalue = fdrtool(rt$pvalue, statistic = "pvalue")$qval   

gene=rownames(rt)
hr=sprintf("%.3f",rt$"HR")
hrLow=sprintf("%.3f",rt$"HR.95L")
hrHigh=sprintf("%.3f",rt$"HR.95H")
Hazard.ratio=paste0(hr,"(",hrLow,"-",hrHigh,")")

pVal=ifelse(rt$qvalue<0.001, "<0.001", sprintf("%.3f", rt$qvalue))

pdf(file=outFile, width = 6, height =6)
n=nrow(rt)
nRow=n+1
ylim=c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2))

xlim = c(0,3)
par(mar=c(4,2,1.5,1.5))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex);text(0.5,n+1,'Cancer',cex=text.cex,font=2,adj=1,)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,expression(bold(italic(q) - value)),cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

par(mar=c(4,1,1.5,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.03,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, '#6495ED', '#FF8C00')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()