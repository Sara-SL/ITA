## CorrelationTissue.R

# Set path ----------------------------------------------------------------

myPath = "~/github/ITA"
myDataPath = "~/github/ITA/Data/"

# Load data ---------------------------------------------------------------

setwd(myPath)
genes.filtered = sort(readLines("genes.filtered.txt"))
bloodData = c("ExprM.GSE102459.2020_10_16.gz", "ExprM.GSE128224.2020_10_22.gz", "ExprM.GSE127792.2020_10_22.gz", "ExprM.GSE97590.2020_10_22.gz", "ExprM.GSE86627.2020_10_16.gz", "ExprM.GSE86331.2020_10_16.gz", "ExprM.GSE83951.2020_10_16.gz")
PBMCData = c("ExprM.GSE82152.2020_10_16.gz", "ExprM.GSE120115.2020_10_22.gz", "ExprM.GSE107990.2020_10_16.gz", "ExprM.GSE87186.2020_10_22.gz", "ExprM.GSE59743.2020_10_16.gz", "ExprM.GSE74816.2020_10_16.gz")


# Create expression & correlation matrix for each dataset -------------------------------

#Select tissue: bloodData or PBMCData
tissue = PBMCData

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
  corrL = append(corrL, list(cor(Mx, use="pairwise.complete.obs")))
}


# Statistics --------------------------------------------------------------
install.packages("bigmemory")
library(bigmemory)

# Initiate matrix for t-values
corrMx <- big.matrix(nrow = 14062, ncol = 14062, type = "double", init = 0, backingfile = "corrMx_backing.bin", descriptorfile = "corrMx_descriptor.desc" )
options(bigmemory.allow.dimnames=TRUE)
colnames(corrMx) = genes.filtered
rownames(corrMx) = genes.filtered

# Create correlation table for each gene and calculate statistics
for(ii in 1:length(genes.filtered)){
  CorrTable.GeneX = matrix(nrow=length(tissue), ncol=length(genes.filtered))
  colnames(CorrTable.GeneX) = genes.filtered

  for(jj in 1:length(corrL)){ 
    if(genes.filtered[ii] %in% colnames(corrL[[jj]])){
      CorrTable.GeneX[jj, which(genes.filtered %in% colnames(corrL[[jj]]))] = corrL[[jj]] [, which(colnames(corrL[[jj]]) == genes.filtered[ii])]
    }
  }
  
  #Calculate statistics for the gene 
  mean.x = colMeans(CorrTable.GeneX, na.rm=T)
  n.x = colSums(!is.na(CorrTable.GeneX),na.rm=T)
  sd.x = apply(FUN=sd, X=CorrTable.GeneX, MARGIN=2, na.rm=T)
  t.vector = mean.x*sqrt(n.x)/sd.x
  
  #Save t-values in matrix
  corrMx[ii,] <- t.vector
  print(ii)
}

write.big.matrix(corrMx, "corrMx_PBMC.csv", row.names = TRUE, col.names = TRUE, sep=',')



# Find interactions that only exist in one of the tissues  ----------------


CorrM1 = CorrM.PBMC - CorrM.Blood
CorrM2 = CorrM.Blood - CorrM.PBMC

