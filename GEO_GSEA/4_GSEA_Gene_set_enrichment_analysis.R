# https://cloud.tencent.com/developer/article/1078364
# https://www.youtube.com/watch?v=Wu_3v36cvfI&list=PLKOVv6BeGBMTg7Q9cLFUUTWR9XQG0Evez&index=8

rm(list = ls())
library(tidyverse)
load(file = 'GSE17708_DEG.Rdata')

library(clusterProfiler)
library(org.Hs.eg.db)
gene = head(rownames(nrDEG),1000)
gene.df <- bitr(gene, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
head(gene.df)
# https://www.ncbi.nlm.nih.gov/gene/374918

kk <- enrichKEGG(gene         = gene.df$ENTREZID,
           organism     = 'hsa',
           pvalueCutoff = 0.05)
head(kk)[,1:6]
