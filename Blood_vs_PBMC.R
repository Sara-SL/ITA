## CorrelationTissue.R

install.packages("bigmemory")
library(bigmemory)

## Set path ----------------------------------------------------------------

myPath = "~/github/ITA"
myDataPath = "~/github/ITA/Data/"
myBigDataPath = "~/github/ITA/BigData/"

## Load data ---------------------------------------------------------------

setwd(myPath)
genes.filtered = sort(readLines("genes.filtered.txt"))
bloodData = c("ExprM.GSE102459.2020_10_16.gz", "ExprM.GSE128224.2020_10_22.gz", "ExprM.GSE127792.2020_10_22.gz", "ExprM.GSE97590.2020_10_22.gz", "ExprM.GSE86627.2020_10_16.gz", "ExprM.GSE86331.2020_10_16.gz", "ExprM.GSE83951.2020_10_16.gz")
PBMCData = c("ExprM.GSE82152.2020_10_16.gz", "ExprM.GSE120115.2020_10_22.gz", "ExprM.GSE107990.2020_10_16.gz", "ExprM.GSE87186.2020_10_22.gz", "ExprM.GSE59743.2020_10_16.gz", "ExprM.GSE74816.2020_10_16.gz")


## Create expression & correlation matrix for each dataset -------------------------------

# Blood -------------------------------------------------------------------
setwd(myDataPath)
corrL.blood = list()
for(i in 1:length(bloodData)){
  print(i)
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
  #Only keep the filtered genes
  Mx <- Mx[, -remove_3]
  
  # Make a correlation matrix from the expression matrix
  corrL.blood = append(corrL.blood, list(cor(Mx, use="pairwise.complete.obs")))
}

setwd(myBigDataPath)
# Initiate matrix for t-values
corrMx.blood <- big.matrix(nrow = 14062, ncol = 14062, type = "double", init = NULL, dimnames = genes.filtered, backingfile = "corrMx_blood_backing.bin", descriptorfile = "corrMx_blood_descriptor.desc", shared = TRUE  )

# Create correlation table for each gene and calculate statistics
for(ii in 1:length(genes.filtered)){
  print(ii)
  CorrTable.GeneX = matrix(nrow=length(bloodData), ncol=length(genes.filtered))
  colnames(CorrTable.GeneX) = genes.filtered
  
  for(jj in 1:length(corrL.blood)){ 
    if(genes.filtered[ii] %in% colnames(corrL.blood[[jj]])){
      CorrTable.GeneX[jj, which(genes.filtered %in% colnames(corrL.blood[[jj]]))] = corrL.blood[[jj]] [, which(colnames(corrL.blood[[jj]]) == genes.filtered[ii])]
    }
  }
  
  #Calculate statistics for the gene 
  mean.x = colMeans(CorrTable.GeneX, na.rm=T)
  n.x = colSums(!is.na(CorrTable.GeneX),na.rm=T)
  sd.x = apply(FUN=sd, X=CorrTable.GeneX, MARGIN=2, na.rm=T)
  t.vector = mean.x*sqrt(n.x)/sd.x
  
  #Save t-values in matrix
  corrMx.blood[ii,] <- t.vector
}

# remove unessecary data to save memory 
rm(corrL.blood, L1, L2, Mx, remove_3, i, j, ii, jj, CorrTable.GeneX, mean.x, n.x, sd.x, t.vector)


# PBMC --------------------------------------------------------------------
setwd(myDataPath)
corrL.PBMC = list()
for(i in 1:length(PBMCData)){
  print(i)
  L1 = readLines(gzfile(PBMCData[i])); closeAllConnections()
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
  corrL.PBMC = append(corrL.PBMC, list(cor(Mx, use="pairwise.complete.obs")))
}

setwd(myBigDataPath)
corrMx.PBMC <- big.matrix(nrow = 14062, ncol = 14062, type = "double", init = NULL, dimnames = genes.filtered, backingfile = "corrMx_PBMC_backing.bin", descriptorfile = "corrMx_PBMC_descriptor.desc", shared = TRUE )

# Create correlation table for each gene and calculate statistics
for(ii in 1:length(genes.filtered)){
  print(ii)
  CorrTable.GeneX = matrix(nrow=length(PBMCData), ncol=length(genes.filtered))
  colnames(CorrTable.GeneX) = genes.filtered
  
  for(jj in 1:length(corrL.PBMC)){ 
    if(genes.filtered[ii] %in% colnames(corrL.PBMC[[jj]])){
      CorrTable.GeneX[jj, which(genes.filtered %in% colnames(corrL.PBMC[[jj]]))] = corrL.PBMC[[jj]] [, which(colnames(corrL.PBMC[[jj]]) == genes.filtered[ii])]
    }
  }
  
  #Calculate statistics for the gene 
  mean.x = colMeans(CorrTable.GeneX, na.rm=T)
  n.x = colSums(!is.na(CorrTable.GeneX),na.rm=T)
  sd.x = apply(FUN=sd, X=CorrTable.GeneX, MARGIN=2, na.rm=T)
  t.vector = mean.x*sqrt(n.x)/sd.x
  
  #Save t-values in matrix
  corrMx.PBMC[ii,] <- t.vector
}

# remove unessecary data to save memory  
rm(corrL.PBMC, L1, L2, Mx, remove_3, i, j, ii, jj, CorrTable.GeneX, mean.x, n.x, sd.x, t.vector)

### Statistics --------------------------------------------------------------

## Find interactions that only exist in one of the tissues  ----------------

# Only in blood -----------------------------------------------------------

CorrM1 = matrix(nrow = length(genes.filtered), ncol = length(genes.filtered))
rownames(CorrM1) = genes.filtered
colnames(CorrM1) = genes.filtered
for(i in 1:length(genes.filtered)){
  CorrM1[i,] = corrMx.PBMC[i,]-corrMx.blood[i,]
}

N = 3 #first N interactions
CorrL1 = list()
for(i in 1:ncol(CorrM1)){
  vx = CorrM1[-i,i]
  vx2 = sort(vx, decreasing=T)
  CorrL1[[i]] = rbind(rep(rownames(CorrM1)[i],N),names(vx2[1:N]),vx2[1:N])
}
TopCorrM1 = matrix(unlist(CorrL1),ncol=3,byrow=T)


# Only in PBMC ------------------------------------------------------------

CorrM2 = matrix(nrow = length(genes.filtered), ncol = length(genes.filtered))
rownames(CorrM2) = genes.filtered
colnames(CorrM2) = genes.filtered
for(i in 1:length(genes.filtered)){
  CorrM2[i,] = corrMx.PBMC[i,]-corrMx.blood[i,]
}

CorrL2 = list()
for(i in 1:ncol(CorrM2)){
  vx = CorrM2[-i,i]
  vx2 = sort(vx, decreasing=T)
  CorrL2[[i]] = rbind(rep(rownames(CorrM2)[i],N),names(vx2[1:N]),vx2[1:N])
}
TopCorrM2 = matrix(unlist(CorrL2),ncol=3,byrow=T)


## Create graph ------------------------------------------------------------


install.packages("igraph")
library(igraph)
Corr.Graph1 = simplify(graph.data.frame(TopCorrM1, directed=F), edge.attr.comb="max")
Corr.Graph2 = simplify(graph.data.frame(TopCorrM2, directed=F), edge.attr.comb="max")



## Analysis ----------------------------------------------------------------

# What kind of genes are hub-genes? i.e genes with many interactions
top_genes = sort(degree(Corr.Graph1),decreasing=T)[1:20]
top_genes = sort(degree(Corr.Graph2),decreasing=T)[1:20]

# Measure of centrality based on shortest path
# how much a gene is affected by its surronding nodes
betweeness_genes = betweeness(Corr.Graph1, v = V(Corr.Graph1), directed = FALSE)
betweeness_genes = betweeness(Corr.Graph2, v = V(Corr.Graph2), directed = FALSE)

# Measures how many steps is required to access every other vertex from a given vertex.
# How independent a genes if from its surrounding genes
closenesse_genes = closenesse(Corr.Graph1, vids = V(Corr.Graph))
closenesse_genes = closenesse(Corr.Graph2, vids = V(Corr.Graph))






