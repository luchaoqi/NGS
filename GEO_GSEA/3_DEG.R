# http://bio-info-trainee.com/bioconductor_China/software/limma.html

rm(list = ls())
suppressPackageStartupMessages(library(CLL))
data(sCLLex)
exprSet=exprs(sCLLex)   ##sCLLex是依赖于CLL这个package的一个对象
samples=sampleNames(sCLLex)
pdata=pData(sCLLex)
# *
# get group_list
group_list=as.character(pdata[,2])
dim(exprSet)

# QC
if(F){
  par(cex = 0.7)
  n.sample=ncol(exprSet)
  if(n.sample>40) par(cex = 0.5)
  cols <- rainbow(n.sample*1.2)
  boxplot(exprSet, col = cols,main="expression value",las=2)
}


# *
t = group_list
design <- model.matrix(~0+factor(t))
colnames(design)=levels(factor(t))
rownames(design)=colnames(exprSet)
design


contrast.matrix<-makeContrasts(paste0(unique(t),collapse = "-"),levels = design)
contrast.matrix ##这个矩阵声明，我们要把progres.组跟stable进行差异分析比较


# limma

##step1
fit <- lmFit(exprSet,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要，大家可以自行看看效果
fit2 <- eBayes(fit2)  ## default no trend !!!
##eBayes() with trend=TRUE
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(nrDEG)


###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

# https://github.com/jmzeng1314/GEO/blob/565eeb9e25f4aee4c5f77f13338cbfc525e1fdc2/task1-check-specific-genes/step3-DEG.R


rm(list = ls())
load(file='GSE17708_new_exprSet.Rdata')
exprSet=new_exprSet
dim(exprSet)
group_list

# select certain genes
exprSet=exprSet[,c(1:3,24:26)]
group_list=group_list[c(1:3,24:26)]

library(limma)
if(F){
  tmp=data.frame(case=c(0,0,0,1,1,1),
               control=c(1,1,1,0,0,0))
}
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)
design

contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),
                               levels = design)
contrast.matrix<-makeContrasts("group72h-grouph",
                               levels = design)

contrast.matrix 
##这个矩阵声明，我们要把72h处理组跟未处理的细胞系进行差异分析比较

deg = function(exprSet,design,contrast.matrix){
  ##step1
  fit <- lmFit(exprSet,design)
  ##step2
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  ##这一步很重要，大家可以自行看看效果
  
  fit2 <- eBayes(fit2)  ## default no trend !!!
  ##eBayes() with trend=TRUE
  ##step3
  tempOutput = topTable(fit2, coef=1, n=Inf)
  nrDEG = na.omit(tempOutput) 
  #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
  head(nrDEG)
  return(nrDEG)
}

re = deg(exprSet,design,contrast.matrix)
nrDEG=re

## heatmap
library(pheatmap)
choose_gene=head(rownames(nrDEG),50) ## 50 maybe better
choose_matrix=exprSet[choose_gene,]
choose_matrix=t(scale(t(choose_matrix)))
pheatmap(choose_matrix,filename = 'DEG_top50_heatmap.png')


library(ggplot2)
## volcano plot
colnames(nrDEG)
plot(nrDEG$logFC,-log10(nrDEG$P.Value))

DEG=nrDEG


logFC_cutoff <- with(DEG,mean(abs( logFC)) + 2*sd(abs( logFC)) )
# logFC_cutoff=1

DEG$change = as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(DEG[DEG$change =='UP',]) ,
                    '\nThe number of down gene is ',nrow(DEG[DEG$change =='DOWN',])
)

g = ggplot(data=DEG, 
           aes(x=logFC, y=-log10(P.Value), 
               color=change)) +
  geom_point(alpha=0.4, size=1.75) +
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle( this_tile ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red')) ## corresponding to the levels(res$change)
print(g)
ggsave(g,filename = 'volcano.png')

save(new_exprSet,group_list,nrDEG,DEG, 
     file='GSE17708_DEG.Rdata')
