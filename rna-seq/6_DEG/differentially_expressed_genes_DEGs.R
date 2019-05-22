# Thanks to jimmy's help
# multiple treat groups:
# https://raw.githubusercontent.com/jmzeng1314/my-R/master/10-RNA-seq-3-groups/hisat2_mm10_htseq.R
# https://github.com/jmzeng1314/GEO/blob/master/task4-NPC/step3-DEG.R
# treat vs. untreat


# BiocManager::install(c("airway","DESeq2","edgeR","limma"))
#What this means is that if you supply a character string to library(),
#you must set character.only to TRUE.

Packages = c("airway","DESeq2","edgeR","limma")
lapply(Packages, library, character.only = TRUE)
rm(Packages)



if(T){
  library(airway)
  data(airway)
  exprSet = assay(airway)
  group_list = colData(airway)[,3]
  save(exprSet, file = 'airway_exprSet.Rdata')
}

# rm(list = ls())
options(stringsAsFactors = F)
load(file = 'airway_exprSet.Rdata')



######################################################################
###################      Firstly for DEseq2      #####################
######################################################################
suppressMessages(library(DESeq2))
colData = data.frame(row.names = colnames(exprSet),
                     group_list = group_list)
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = colData,
                              design = ~ group_list)
dds <- DESeq(dds)

# > group_list
# [1] untrt trt   untrt trt   untrt trt   untrt trt  

res <- results(dds,
               contrast=c("group_list","trt","untrt"))
resOrdered <- res[order(res$padj),]
head(resOrdered)
DEG =as.data.frame(resOrdered)
write.csv(DEG,"DEG.results.csv")



DEG = na.omit(DEG)
if(F){
  png("qc_dispersions.png", 1000, 1000, pointsize=20)
  plotDispEsts(dds, main="Dispersion plot")
  dev.off()
  
  rld <- rlogTransformation(dds)
  exprMatrix_rlog=assay(rld) 
  write.csv(exprMatrix_rlog,'exprMatrix.rlog.csv' )
  
  
  
  png("DEseq_RAWvsNORM.png",height = 800,width = 800)
  par(cex = 0.7)
  n.sample=ncol(exprSet)
  if(n.sample>40) par(cex = 0.5)
  cols <- rainbow(n.sample*1.2)
  par(mfrow=c(2,2))
  boxplot(exprSet, col = cols,main="expression value",las=2)
  boxplot(exprMatrix_rlog, col = cols,main="expression value",las=2)
  hist(as.matrix(exprSet))
  hist(exprMatrix_rlog)
  dev.off()
  
}



nrDEG = DEG
DEGseq_DEG = DEG

# heatmap/volcano
if(F){
  library(pheatmap)
  choose_gene=head(rownames(nrDEG),50) ## 50 maybe better
  choose_matrix = exprSet[choose_gene,]
  choose_matrix=t(scale(t(choose_matrix)))
  pheatmap(choose_matrix,filename = 'DEG_top50_heatmap.png')
  
  
  
  # volcano
  
  logFC_cutoff <- with(DEG,mean(abs( log2FoldChange)) + 2*sd(abs( log2FoldChange)) )
  # logFC_cutoff=1
  
  DEG$change = as.factor(ifelse(DEG$pvalue < 0.05 & abs(DEG$log2FoldChange) > logFC_cutoff,
                                ifelse(DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT')
  )
  this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                      "\nThe number of up gene is ",nrow(DEG[DEG$change =='UP',]) ,
                      "\nThe number of down gene is ",nrow(DEG[DEG$change =='DOWN',])
  )
  
  library(ggplot2)
  g = ggplot(data=DEG, 
             aes(x=log2FoldChange, y=-log10(pvalue), 
                 color=change)) +
    geom_point(alpha=0.4, size=1.75) +
    theme_set(theme_set(theme_bw(base_size=20)))+
    xlab("log2 fold change") + ylab("-log10 p-value") +
    ggtitle( this_tile ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
    scale_colour_manual(values = c('blue','black','red')) ## corresponding to the levels(res$change)
  print(g)
  ggsave(g,filename = 'volcano.png')
}


######################################################################
###################      Then  for edgeR        #####################
######################################################################

if(T){
  
  library(edgeR)
  d <- DGEList(counts=exprSet,group=factor(group_list))
  d$samples$lib.size <- colSums(d$counts)
  d <- calcNormFactors(d)
  d$samples
  
  ## The calcNormFactors function normalizes for RNA composition by finding a set of scaling
  ## factors for the library sizes that minimize the log-fold changes between the samples for most
  ## genes. The default method for computing these scale factors uses a trimmed mean of Mvalues
  ## (TMM) between each pair of samples
  
  png('edgeR_MDS.png')
  plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
  dev.off()
  
  # The glm approach to multiple groups is similar to the classic approach, but permits more general comparisons to be made
  
  dge=d
  
  design <- model.matrix(~0+factor(group_list))
  rownames(design)<-colnames(dge)
  colnames(design)<-levels(factor(group_list))
  
  dge <- estimateGLMCommonDisp(dge,design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  
  fit <- glmFit(dge, design)
  
  
  # > design
  # trt untrt
  # SRR1039508   0     1
  # SRR1039509   1     0
  # SRR1039512   0     1
  # SRR1039513   1     0
  # SRR1039516   0     1
  # SRR1039517   1     0
  # SRR1039520   0     1
  # SRR1039521   1     0
  

  lrt <- glmLRT(fit,  contrast=c(1,0))
  nrDEG=topTags(lrt, n=nrow(exprSet))
  nrDEG=as.data.frame(nrDEG)
  head(nrDEG)
  edgeR_DEG = nrDEG
  # write.csv(nrDEG,"DEG_treat_12_edgeR.csv",quote = F)
  # 
  # lrt <- glmLRT(fit, contrast=c(0,1) )
  # nrDEG=topTags(lrt, n=nrow(exprSet))
  # nrDEG=as.data.frame(nrDEG)
  # head(nrDEG)
  # write.csv(nrDEG,"DEG_treat_2_edgeR.csv",quote = F)
  # summary(decideTests(lrt))
  # plotMD(lrt)
  # abline(h=c(-1, 1), col="blue")
}

######################################################################
###################      Then  for limma/voom        #################
######################################################################

if(T){
  suppressMessages(library(limma))
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(exprSet)
  
  dge <- DGEList(counts=exprSet)
  dge <- calcNormFactors(dge)
  logCPM <- cpm(dge, log=TRUE, prior.count=3)
  
  v <- voom(dge,design,plot=TRUE, normalize="quantile")
  fit <- lmFit(v, design)
  
  group_list
  # cont.matrix=makeContrasts(contrasts=c('treat_12-control','treat_2-control'),levels = design)
  cont.matrix=makeContrasts(contrasts=c('untrt-trt'),levels = design)
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)
  tempOutput = topTable(fit2, coef='untrt-trt', n=Inf)
  DEG_limma_voom = na.omit(tempOutput)
  save(DEG_limma_voom, DEGseq_DEG, edgeR_DEG, dds, exprSet, group_list, file = 'DEG_results.Rdata')
  
  if(F){
    write.csv(DEG_limma_voom,"DEG_limma_voom.csv",quote = F)
    png("limma_voom_RAWvsNORM.png",height = 600,width = 600)
    exprSet_new=v$E
    par(cex = 0.7)
    n.sample=ncol(exprSet)
    if(n.sample>40) par(cex = 0.5)
    cols <- rainbow(n.sample*1.2)
    par(mfrow=c(2,2))
    boxplot(exprSet, col = cols,main="expression value",las=2)
    boxplot(exprSet_new, col = cols,main="expression value",las=2)
    hist(as.matrix(exprSet))
    hist(exprSet_new)
    dev.off()
  }
}



######################################################################
##################      Futher investigation       ###################
######################################################################

# Compare all the results from DEseq2, edgeR and limma

rm(list = ls())
load(file = 'DEG_results.Rdata')


a1 = data.frame(gene = rownames(DEG_limma_voom),
                limma = DEG_limma_voom$logFC)

a2 = data.frame(gene = rownames(DEGseq_DEG),
                DEseq = DEGseq_DEG$log2FoldChange)

a3 = data.frame(gene = rownames(edgeR_DEG),
                edgeR = edgeR_DEG$logFC)



