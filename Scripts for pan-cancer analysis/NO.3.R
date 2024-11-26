library(plyr)
library(reshape2)
library(ggpubr)
inputFile=" UCHL5geneExp.txt"      
setwd("pathway")    

data=read.table(inputFile, header=T, sep="\t", check.names=F)
gene=colnames(data)[2]
colnames(data)[2]="expression"

p=ggboxplot(data, x="CancerType", y="expression", color="Type",
     xlab="",
     ylab=paste0(gene," expression"),
     width=0.6,
     palette = c("#6495ED","#FF8C00") )
p=p+rotate_x_text(50)

tumor_normal = data.frame(table(unique(data[,c('Type','CancerType')])$CancerType))
tumor_normal = tumor_normal[which(tumor_normal$Freq == 2),]

data2 = data[data$CancerType %in% tumor_normal$Var1,]

temp = data.frame()
for (i in 1:length(unique(data2$CancerType))) {
    dat = data2[which(data2$CancerType == unique(data2$CancerType)[i]),]
	normal = dat[which(dat$Type == 'Normal'),]
	tumor = dat[which(dat$Type == 'Tumor'),]
	wilcox_results = wilcox.test(tumor$expression, normal$expression)
	temp[i,1] = unique(data2$CancerType)[i]
	temp[i,2] = wilcox_results$p.value
}
colnames(temp) = c('CancerType','pvalue')
temp$fdr <- p.adjust(temp$pvalue, method = "fdr")
temp = temp[order(temp$fdr),]

cut005 = temp[which(temp$fdr <= 0.05 & temp$fdr > 0.01),]
dim(cut005)[2] 
cut001 = temp[which(temp$fdr <= 0.01 & temp$fdr >= 0.001),]
dim(cut005)[1]
cut0001 = temp[which(temp$fdr <= 0.001),]
dim(cut0001)[1]

cut0001 = max(cut0001$pvalue) + 10^-20
cut001 = max(cut001$pvalue) + 10^-20
cut005 = 0.05

p=p+stat_compare_means(aes(group=Type),
      method="wilcox.test",
      symnum.args=list(cutpoints = c(0, cut0001, cut001, cut005, 1), symbols = c("***", "**", "*", " ")),
      label = "p.signif")
pdf(file="diff.pdf", width=8, height=5.5)
print(p)
dev.off()
