library(plyr)
library(reshape2)
library(ggpubr)

inputFile = "UCHL5geneExp.txt"     
setwd("pathway")    

data = read.table(inputFile, header = T, sep = "\t", check.names = F)
gene = colnames(data)[2]
colnames(data)[2] = "expression"

data = data[(data[, "Type"] == "Tumor"),]
med = ddply(data, "CancerType", summarise, med = median(expression))
data$CancerType = factor(data$CancerType, levels = med[order(med[, "med"], decreasing = T), "CancerType"])

# **Kruskal-Wallis test**
kruskal_result = kruskal.test(expression ~ CancerType, data = data)
print(kruskal_result)

p = ggboxplot(data, 
              x = "CancerType", 
              y = "expression", 
              fill = "CancerType",
              xlab = "",
              ylab = paste0(gene, " expression"),
              width = 0.6,
              legend = "none") + 
    rotate_x_text(50)

p = p + labs(subtitle = paste0("Kruskal-Wallis p-value: ", signif(kruskal_result$p.value, 3)))

pdf(file = "boxplot_with_kruskal.pdf", width = 8, height = 5.5)
print(p)
dev.off().