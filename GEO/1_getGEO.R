# https://github.com/jmzeng1314/GEO

# eSet[["GSE3325_series_matrix.txt.gz"]]@annotation
# This is GPL570 platform

# google it or find package here: http://www.bio-info-trainee.com/1399.html
# https://github.com/jmzeng1314/my-R/blob/master/1-get-all-probeset-info/GPL_info.csv
# then we use package library(hgu133plus2.db) in step2_id_filter

rm(list=ls())
if(F){
  library(GEOquery)
  eSet <- getGEO('GSE3325', destdir=".",
                 AnnotGPL = F,
                 getGPL = F)
  save(eSet,file='GSE3325_eSet.Rdata')
}
load('GSE3325_eSet.Rdata')
b = eSet[[1]]

# raw_exprSet = eSet[["GSE3325_series_matrix.txt.gz"]]@assayData[["exprs"]]
raw_exprSet=exprs(b) 
phe=pData(b)
library(stringr)
group_list= str_split(as.character(phe$title),' ',simplify = T)[,1]
save(raw_exprSet,group_list,
     file='GSE3325_raw_exprSet.Rdata')



rm(list=ls())
if(F){
  library(GEOquery)
  eSet <- getGEO('GSE17708', destdir=".",
                 AnnotGPL = F,
                 getGPL = F)
  save(eSet,file='GSE17708_eSet.Rdata')
}
load('GSE17708_eSet.Rdata')
b = eSet[[1]]
raw_exprSet=exprs(b) 
raw_exprSet[1:4,1:4]
phe=pData(b)
phe$title
library(stringr)
group_list= str_split(as.character(phe$title),' ',simplify = T)[,11]
group_list=paste0('group',group_list,'h')
save(raw_exprSet,group_list,
     file='GSE17708_raw_exprSet.Rdata')
