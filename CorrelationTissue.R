## CorrelationTissue.R

# Set path ----------------------------------------------------------------

myPath = "/Users/aina/Documents/X/Slutkurs/github/ITA"
myDataPath = "/Users/aina/Documents/X/Slutkurs/github/ITA/Data/"

# Load data ---------------------------------------------------------------

setwd(myPath)
<<<<<<< HEAD
genes.filtered = sort(readLines("/Users/aina/Documents/X/Slutkurs/github/ITA/genes.filtered.txt"))
=======
genes.filtered = sort(readLines("genes.filtered.txt"))
>>>>>>> 94758c2aef536b0728d3bddd6f938656518f1b8b
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

# Progress bar -----------------------------------------------------------
install.packages("svMisc")
library(svMisc)
require(svMisc)

# Statistics --------------------------------------------------------------
<<<<<<< HEAD
# Extract correlation values for each gene
corrMx = matrix(nrow = 0, ncol=length(genes.filtered))
colnames(corrMx) = genes.filtered
install.packages("bigmemory")
library(bigmemory)
corrMx <- big.matrix(nrow = 14062, ncol = 14062, type = "double", init = 0, backingfile = "corrMx_backing.bin", descriptorfile = "corrMx_descriptor.desc" )
options(bigmemory.allow.dimnames=TRUE)
colnames(corrMx) = genes.filtered
rownames(corrMx) = genes.filtered

# Create correlation table for the gene and calculate statistics
for(ii in 1:length(genes.filtered)){
  CorrTable.GeneX = matrix(nrow=length(tissue), ncol=length(genes.filtered))
  colnames(CorrTable.GeneX) = genes.filtered
  
  for(jj in 1:length(corrL)){ 				    	 
    if(genes.filtered[ii] %in% colnames(corrL[[jj]])){
      CorrTable.GeneX[jj, which(genes.filtered %in% colnames(corrL[[jj]]))] = corrL[[jj]] [, which(colnames(corrL[[jj]]) == genes.filtered[ii])]
    } else {
      next
    }
  }
        
  #Calculate statistics for the gene 
  mean.x = colMeans(CorrTable.GeneX, na.rm=TRUE)
  n.x = colSums(!is.na(CorrTable.GeneX),na.rm=TRUE)
  sd.x = apply(FUN=sd, X=CorrTable.GeneX, MARGIN=2, na.rm=TRUE)
  t.vector = mean.x*sqrt(n.x)/sd.x
        
  corrMx[ii,] <- t.vector
  print(ii)
}
corrMx_blood <- corrMx

# Read corrMx for PBMC
setwd('/Users/aina/Documents/X/Slutkurs')
corrMx_PBMC = read.big.matrix("corrMx_PBMC.csv", sep = ",", header = TRUE, has.row.names = TRUE, descriptorfile = "corrMx_descriptor.desc", backingfile = "corrMx_backing.bin")

# Find interactions that only exist in one of the tissues  --------------------------------------------------------------
# CorrM1 repsresents PBMC-blood 
CorrM1 = matrix(nrow = length(genes.filtered), ncol = length(genes.filtered))
rownames(CorrM1) = genes.filtered
colnames(CorrM1) = genes.filtered

for(i in 1:length(genes.filtered)){
  CorrM1[i,] = corrMx_PBMC[i,]-corrMx_blood[i,]
}

# Corr M2 represents blood-PBMC
CorrM2 = matrix(nrow = length(genes.filtered), ncol = length(genes.filtered))
rownames(CorrM2) = genes.filtered
colnames(CorrM2) = genes.filtered

for(i in 1:length(genes.filtered)){
  CorrM2[i,] = corrMx_blood[i,]-corrMx_PBMC[i,]
}

#Select matrix: CorrM1, CorrM2 ------
CorrM = CorrM1

N = 3 #first N interactions
CorrL = list()
for(i in 1:ncol(CorrM)){
  vx = CorrM[-i,i]
  vx2 = sort(vx, decreasing=T)
  CorrL[[i]] = rbind(rep(rownames(CorrM)[i],N),names(vx2[1:N]),vx2[1:N])
}
TopCorrM = matrix(unlist(CorrL),ncol=3,byrow=T)
colnames(TopCorrM) = c("Gene 1", "Gene 2", "T-value")

hist(as.numeric(TopCorrM[,3]), xlim = c(0,100), breaks = length(TopCorrM[,3]))

# Cutoff = 50
remove = c()
for(j in 1:length(TopCorrM[,3])){
  if(as.numeric(TopCorrM[j,3]) < 50){
    remove = append(remove, j)
  } else{
    next
  }
}

TopCorrM_filtered = TopCorrM[-remove,]

# Create graph -----
install.packages("igraph")
library(igraph)

Corr.Graph = simplify(graph.data.frame(TopCorrM_filtered[,1:2], directed=FALSE), edge.attr.comb="max")
plot.igraph(Corr.Graph, vertex.size=3, vertex.label=NA)

# Analysis ----
# What kind of genes are hub-genes? i.e genes with many interactions
top_genes = sort(degree(Corr.Graph),decreasing=TRUE)[1:20]

# Plot hug-gene interactions
HubGenes = which(TopCorrM_filtered[,1:2] %in% names(top_genes))
Corr.Graph_2 = simplify(graph.data.frame(TopCorrM_filtered[HubGenes[1:200], 1:2], directed=FALSE), edge.attr.comb="max")
plot.igraph(Corr.Graph_2)

# Measure of centrality based on shortest path
# how much a gene is affected by its surronding nodes
betweenness_genes = betweenness(Corr.Graph, v = V(Corr.Graph), directed = FALSE)

# Measures how many steps is required to access every other vertex from a given vertex.
# How independent a genes if from its surrounding genes
closeness_genes = closeness(Corr.Graph, vids = V(Corr.Graph))
=======
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

>>>>>>> 94758c2aef536b0728d3bddd6f938656518f1b8b
