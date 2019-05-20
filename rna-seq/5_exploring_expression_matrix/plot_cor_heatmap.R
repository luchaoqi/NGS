rm(list = ls())
options(stringsAsFactors = F)

a = read.table('all.id.txt',header = T)
meta = a[,1:6]
exprSet = a[,7:ncol(a)]
a1 = exprSet[,'SRR1039516']
colnames(exprSet)
group_list = colData(airway)[,3]


library(corrplot)
png(cor.png)
corrplot(cor(log2(exprSet+1)))
dev.off()



library(pheatmap)
png(heatmap.png)
m = cor(log2(exprSet+1))
pheatmap(scale(m))
dev.off()



a1 = data.frame(id = names(a1),a1=as.numeric(a1))
a2 = data.frame(id = meta[,1],a2 = a2)
a2$id = strsplit(a2$id,'\\.',simplify=T)[,1]
temp = merge(a1,a2,by = 'id')
plot(temp[,2:3])





# http://bioconductor.org/install/

# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install()

library(airway)
data(airway)
exprSet = assay(airway)
group_list = colData(airway)[,3]

## hclust
colnames(exprSet)=paste(group_list,1:ncol(exprSet),sep='_')
# Define nodePar
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19),
                cex = 0.7, col = "blue")
hc=hclust(dist(t(log2(exprSet+1))))
par(mar=c(5,5,5,10))
png('hclust.png',res=120)
plot(as.dendrogram(hc), nodePar = nodePar, horiz = TRUE)
dev.off()