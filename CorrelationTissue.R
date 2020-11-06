## CorrelationTissue.R

# Set path ----------------------------------------------------------------

myPath = "~/github/ITA"
myDataPath = "~/github/ITA/Data/"

# Load data ---------------------------------------------------------------

setwd(myPath)
genes.filtered = sort(readLines("~/github/ITA/genes.filtered.txt"))
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
  
  Mx <- Mx[, -remove_3]
  
  # Make a correlation matrix from the expression matrix
  corrL = append(corrL, list(cor(Mx, use="pairwise.complete.obs")))
}


# Statistics --------------------------------------------------------------

CorrTable.GeneX = matrix(nrow=X1, ncol=X2)
colnames(CorrTable.GeneX) = genes.filtered

for(ii in 1:length(genes.filtered)){
  for(jj in 1:length(CorrL_PBMC)){
    if(genes.filtered[ii]%in%colnames(CorrL_PBMC[[jj]])){
      Corr.Table.GeneX[jj,] = CorrL_PBMC[[jj]][,genes.filtered[ii]]
    }
  }
}
