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