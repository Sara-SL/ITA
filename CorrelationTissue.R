## Blood dataset
# Create expression matrix for each dataset
setwd("~/Documents/X/Slutkurs/github/ITA")

genes.filtered = sort(readLines("~/Documents/X/Slutkurs/github/ITA/genes.filtered.txt"))
bloodData = c("ExprM.GSE102459.2020_10_16.gz", "ExprM.GSE128224.2020_10_22.gz", "ExprM.GSE127792.2020_10_22.gz", "ExprM.GSE97590.2020_10_22.gz", "ExprM.GSE86627.2020_10_16.gz", "ExprM.GSE86331.2020_10_16.gz", "ExprM.GSE83951.2020_10_16.gz")

setwd("~/Documents/X/Slutkurs/github/ITA/Data")
corrL = list()
for(i in 1:length(bloodData)){
  L1 = readLines(gzfile(bloodData[i])); closeAllConnections()
  L2 = strsplit(L1, split=",")
  Mx = sapply(FUN=function(v){return(as.numeric(v[-1]))},X=L2[2:length(L2)])
  rownames(Mx) = L2[[1]][-1]
  colnames(Mx) = sapply(FUN=getElement, X=L2[2:length(L2)], 1)
  for(i in 1:nrow(Mx)){Mx[i,] = Mx[i,]/max(Mx[i,],na.rm=T)}
  
  remove_3 = c()
  for(j in 1:ncol(Mx)){
    if(colnames(Mx)[j] %in% genes.filtered){
      next()
    }else{
      remove_3 = append(remove_3, j)
    }
  }
  
  Mx <- Mx[, -remove_3]
  
  #setwd("~/Documents/X/Slutkurs/github/ITA/Filtered_data/Blood")
  #file_name = paste("fitered_", strsplit(bloodData[i], split=".gz"),".csv", sep='')
  #write.csv(VarMatrix, file = file_name, row.names = TRUE)
  
  # Make a correlation matrixfrom the expression matrix
  corrL = append(corrL, list(cor(Mx, use="pairwise.complete.obs")))
}
