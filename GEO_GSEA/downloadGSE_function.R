# function 
downGSE <- function(studyID = "GSE1009", destdir = ".") {
  
  library(GEOquery)
  eSet <- getGEO(studyID, destdir = destdir, getGPL = F)
  
  exprSet = exprs(eSet[[1]])
  pdata = pData(eSet[[1]])
  
  write.csv(exprSet, paste0(studyID, "_exprSet.csv"))
  write.csv(pdata, paste0(studyID, "_metadata.csv"))
  return(eSet)
  
}
