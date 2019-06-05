
rm(list=ls()) 
load(file='GSE17708_raw_exprSet.Rdata')

exprSet=raw_exprSet


# lib the package based on the platform check step1 for further info
# eSet[["GSE17708_series_matrix.txt.gz"]]@annotation

# google it or find package here: http://www.bio-info-trainee.com/1399.html
# https://github.com/jmzeng1314/my-R/blob/master/1-get-all-probeset-info/GPL_info.csv
# then we use package library(hgu133plus2.db) in step2_id_filter


library(hgu133plus2.db)
ids=toTable(hgu133plus2SYMBOL)
length(unique(ids$symbol))
tail(sort(table(ids$symbol)))
# hist of symbol
table(sort(table(ids$symbol)))
plot(table(sort(table(ids$symbol))))




table(rownames(exprSet) %in% ids$probe_id)
dim(exprSet)
exprSet=exprSet[rownames(exprSet) %in% ids$probe_id,]
dim(exprSet)

ids=ids[match(rownames(exprSet),ids$probe_id),]
head(ids)
exprSet[1:5,1:5]

jimmy <- function(exprSet,ids){
  tmp = by(exprSet,
           ids$symbol,
           function(x) rownames(x)[which.max(rowMeans(x))] )
  probes = as.character(tmp)
  print(dim(exprSet))
  exprSet=exprSet[rownames(exprSet) %in% probes ,]
  
  print(dim(exprSet))
  rownames(exprSet)=ids[match(rownames(exprSet),ids$probe_id),2]
  return(exprSet)
}

new_exprSet <- jimmy(exprSet,ids)

save(new_exprSet,group_list,
     file='GSE42872_new_exprSet.Rdata')


load(file='GSE42872_new_exprSet.Rdata')
exprSet=new_exprSet

# 看下一些常见基因的表达量是不是正常的(高于或者低于一般的表达量)
exprSet['GAPDH',]
boxplot(exprSet[,1])
exprSet['ACTB',]
# 发现高于说明看起来是正常的 之前的id转化应该对的



# 第二个探索就是PCA和hclust下面的代码是基于药物使用后的时间的 可以把grouplist改为untreated和treated
# group_list:
# b = eSet[[1]]
## raw_exprSet = eSet[["GSE3325_series_matrix.txt.gz"]]@assayData[["exprs"]]
# raw_exprSet=exprs(b) 
# phe=pData(b)

if(T){
  
  library(reshape2)
  exprSet_L=melt(exprSet)
  colnames(exprSet_L)=c('probe','sample','value')
  exprSet_L$group=rep(group_list,each=nrow(exprSet))
  head(exprSet_L)
  ### ggplot2
  library(ggplot2)
  # 基本上中间均值线?在同一水平 不然的话就有batch effect
  # sol:sv combine
  p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_boxplot()
  print(p)
  p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_violin()
  print(p)
  p=ggplot(exprSet_L,aes(value,fill=group))+geom_histogram(bins = 200)+facet_wrap(~sample, nrow = 4)
  print(p)
  p=ggplot(exprSet_L,aes(value,col=group))+geom_density()+facet_wrap(~sample, nrow = 4)
  print(p)
  p=ggplot(exprSet_L,aes(value,col=group))+geom_density()
  print(p)
  p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_boxplot()
  p=p+stat_summary(fun.y="mean",geom="point",shape=23,size=3,fill="red")
  p=p+theme_set(theme_set(theme_bw(base_size=20)))
  p=p+theme(text=element_text(face='bold'),axis.text.x=element_text(angle=30,hjust=1),axis.title=element_blank())
  print(p) 
  ## hclust
  colnames(exprSet)=paste(group_list,1:ncol(exprSet),sep='_')
  # Define nodePar
  nodePar <- list(lab.cex = 0.6, pch = c(NA, 19),
                  cex = 0.7, col = "blue")
  hc=hclust(dist(t(exprSet)))
  par(mar=c(5,5,5,10))
  png('hclust.png',res=120)
  plot(as.dendrogram(hc), nodePar = nodePar, horiz = TRUE)
  dev.off()
  
  ## PCA
  
  library(ggfortify)
  df=as.data.frame(t(exprSet))
  df$group=group_list
  png('pca.png',res=120)
  autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')+theme_bw()
  dev.off()
}

