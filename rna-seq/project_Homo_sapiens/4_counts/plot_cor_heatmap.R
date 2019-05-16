rm(list = ls())
options(stringsAsFactors = F)

a = read.table('all.id.txt',header = T)
meta = a[,1:6]
exprSet = a[,7:ncol(a)]
colnames(exprSet)



library(corrplot)
png(cor.png)
corrplot(cor(log2(exprSet+1)))
dev.off()



library(pheatmap)
png(heatmap.png)
m = cor(log2(exprSet+1))
pheatmap(scale(m))
dev.off()



