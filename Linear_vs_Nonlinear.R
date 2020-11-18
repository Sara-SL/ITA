## Mutual information vs Person correlation
# Librarys and packages ---------------------------------------------------
library(infotheo)

# Set path ----------------------------------------------------------------

myPath = "/Users/aina/Documents/X/Slutkurs/github/ITA"
myDataPath = "/Users/aina/Documents/X/Slutkurs/github/ITA/Data/"

# Load data ---------------------------------------------------------------

setwd(myPath)
genes.filtered = sort(readLines("/Users/aina/Documents/X/Slutkurs/github/ITA/genes.filtered.txt"))
bloodData = c("ExprM.GSE102459.2020_10_16.gz", "ExprM.GSE128224.2020_10_22.gz", "ExprM.GSE127792.2020_10_22.gz", "ExprM.GSE97590.2020_10_22.gz", "ExprM.GSE86627.2020_10_16.gz", "ExprM.GSE86331.2020_10_16.gz", "ExprM.GSE83951.2020_10_16.gz")
PBMCData = c("ExprM.GSE82152.2020_10_16.gz", "ExprM.GSE120115.2020_10_22.gz", "ExprM.GSE107990.2020_10_16.gz", "ExprM.GSE87186.2020_10_22.gz", "ExprM.GSE59743.2020_10_16.gz", "ExprM.GSE74816.2020_10_16.gz")


# Create expression & correlation matrix for each dataset -------------------------------
#Select tissue: bloodData or PBMCData
tissue = bloodData

setwd(myDataPath)
corrL = list()
for(i in 1:length(tissue)){
  L1 = readLines(gzfile(tissue[i])); closeAllConnections()
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
  #Only keep the filtered genes
  Mx <- Mx[, -remove_3]
  
  # Make a correlation matrix from the expression matrix
  Mx_disc = discretize(Mx, disc = "equalwidth")
  corrL = append(corrL, list(mutinformation(Mx_disc, method = "mm")))
  
  print(i)
}

for(l in 1:length(corrL)){
  corrL[[l]] = replace(corrL[[l]], corrL[[l]]<=0.7, NA)
  print(l)
}
