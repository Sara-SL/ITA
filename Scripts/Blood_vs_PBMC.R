## Blood_vs_PBMC.R

library(infotheo)
library(bigmemory)
library(igraph)

## Set path ----------------------------------------------------------------
myPath = "~/github/ITA/GeneData/"
myDataPath = "~/github/ITA/ExpressionData/"
myBigDataPath = "~/github/ITA/BigData/"
myTablePath = "~/github/ITA/SupplementaryTables/"
myFigurePath = "~/github/ITA/SupplementaryFigures/"

## Load data ---------------------------------------------------------------
setwd(myPath)
genes.filtered = sort(readLines("genes.filtered_blood_vs_PBMC.txt"))
bloodData = c("ExprM.GSE102459.2020_10_16.gz", "ExprM.GSE128224.2020_10_22.gz", "ExprM.GSE127792.2020_10_22.gz", "ExprM.GSE97590.2020_10_22.gz", "ExprM.GSE86627.2020_10_16.gz", "ExprM.GSE86331.2020_10_16.gz", "ExprM.GSE83951.2020_10_16.gz")
PBMCData = c("ExprM.GSE82152.2020_10_16.gz", "ExprM.GSE120115.2020_10_22.gz", "ExprM.GSE107990.2020_10_16.gz", "ExprM.GSE87186.2020_10_22.gz", "ExprM.GSE59743.2020_10_16.gz", "ExprM.GSE74816.2020_10_16.gz")
geneNames = as.matrix(read.table("hsapiens.SYMBOL.txt", sep="\t", header=F))
ann.v = geneNames[,2]; 
names(ann.v) = geneNames[,1]

#### -------------------- Create expression & correlation matrix for each dataset -------------------------------

# Blood -------------------------------------------------------------------

# Create expression matrix
corrL.blood = list()
for(i in 1:length(bloodData)){
  setwd(myDataPath)
  print(i)
  L1 = readLines(gzfile(bloodData[i])); closeAllConnections()
  L2 = strsplit(L1, split=",")
  Mx = sapply(FUN=function(v){return(as.numeric(v[-1]))},X=L2[2:length(L2)])
  rownames(Mx) = L2[[1]][-1]
  colnames(Mx) = sapply(FUN=getElement, X=L2[2:length(L2)], 1)
  for(i in 1:nrow(Mx)){Mx[i,] = Mx[i,]/max(Mx[i,],na.rm=T)}
  
  #Only keep the filtered genes
  genes.remove = c()
  for(j in 1:ncol(Mx)){
    if(colnames(Mx)[j] %in% genes.filtered){
      next()
    }else{
      genes.remove = append(genes.remove, j)
    }
  }
  Mx <- Mx[, -genes.remove]
  
  # Make a correlation matrix from the expression matrix
  corrL.blood = append(corrL.blood, list(cor(Mx, use="pairwise.complete.obs")))
}

# Remove correlation values below 0.7
for (i in 1:length(corrL.blood)) {
  print(i)
  corrL.blood[[i]] = replace(corrL.blood[[i]], corrL.blood[[i]] <= 0.7, NA)
}

# Initiate matrix for t-values
setwd(myBigDataPath)
corrMx.blood <- big.matrix(nrow = length(genes.filtered), ncol = length(genes.filtered), type = "double", init = NULL, dimnames = genes.filtered, backingfile = "corrMx_blood_backing.bin", descriptorfile = "corrMx_blood_descriptor.desc", shared = TRUE  )
options(bigmemory.allow.dimnames=TRUE)
colnames(corrMx.blood) = genes.filtered
rownames(corrMx.blood) = genes.filtered

# Create correlation table for each gene and calculate the t-value
for(ii in 1:length(genes.filtered)){
  print(ii)
  CorrTable.GeneX = matrix(nrow=length(bloodData), ncol=length(genes.filtered))
  colnames(CorrTable.GeneX) = genes.filtered
  
  # Put genes in same order 
  for(jj in 1:length(corrL.blood)){ 
    if(genes.filtered[ii] %in% colnames(corrL.blood[[jj]])){
      CorrTable.GeneX[jj, match(colnames(corrL.blood[[jj]]), genes.filtered)] = corrL.blood[[jj]] [, match(genes.filtered[ii], colnames(corrL.blood[[jj]]))]
    }
  }
  
  # Only include correlations found in more than 4 datasets 
  index = which(colSums(!is.na(CorrTable.GeneX)) < 4)
  CorrTable.GeneX[,index] = NA
  
  #Calculate statistics for the gene 
  mean.x = colMeans(CorrTable.GeneX, na.rm=T)
  n.x = colSums(!is.na(CorrTable.GeneX),na.rm=T)
  sd.x = apply(FUN=sd, X=CorrTable.GeneX, MARGIN=2, na.rm=T)
  t.vector = mean.x*sqrt(n.x)/sd.x
  
  #Save t-values in matrix
  corrMx.blood[ii,] <- t.vector
}

# Save matrix to .csv file
setwd(myBigDataPath)
write.big.matrix(corrMx.blood, "corrMx_blood.csv", row.names = TRUE, col.names = TRUE, sep=',')


# PBMC --------------------------------------------------------------------

# Create expression matrix
corrL.PBMC = list()
for(i in 1:length(PBMCData)){
setwd(myDataPath)
  print(i)
  L1 = readLines(gzfile(PBMCData[i])); closeAllConnections()
  L2 = strsplit(L1, split=",")
  Mx = sapply(FUN=function(v){return(as.numeric(v[-1]))},X=L2[2:length(L2)])
  rownames(Mx) = L2[[1]][-1]
  colnames(Mx) = sapply(FUN=getElement, X=L2[2:length(L2)], 1)
  for(i in 1:nrow(Mx)){Mx[i,] = Mx[i,]/max(Mx[i,],na.rm=T)}

  #Only keep the filtered genes
  genes.remove = c()
  for(j in 1:ncol(Mx)){
    if(colnames(Mx)[j] %in% genes.filtered){
      next()
    }else{
      genes.remove = append(genes.remove, j)
    }
  }
  Mx <- Mx[, -genes.remove]
  
  # Make a correlation matrix from the expression matrix
  corrL.PBMC = append(corrL.PBMC, list(cor(Mx, use="pairwise.complete.obs")))
}

# Remove correlation values below 0.7
for (i in 1:length(corrL.PBMC)) {
  print(i)
  corrL.PBMC[[i]] = replace(corrL.PBMC[[i]], corrL.PBMC[[i]] <= 0.7, NA)
}

# Initiate matrix for t-values
setwd(myBigDataPath)
corrMx.PBMC <- big.matrix(nrow = length(genes.filtered), ncol = length(genes.filtered), type = "double", init = NULL, dimnames = genes.filtered, backingfile = "corrMx_PBMC_backing.bin", descriptorfile = "corrMx_PBMC_descriptor.desc", shared = TRUE )
options(bigmemory.allow.dimnames=TRUE)
colnames(corrMx.PBMC) = genes.filtered
rownames(corrMx.PBMC) = genes.filtered

# Create correlation table for each gene and calculate the t-value
for(ii in 1:length(genes.filtered)){
  print(ii)
  CorrTable.GeneX = matrix(nrow=length(PBMCData), ncol=length(genes.filtered))
  colnames(CorrTable.GeneX) = genes.filtered
  
  # Put genes in same order 
  for(jj in 1:length(corrL.PBMC)){ 
    if(genes.filtered[ii] %in% colnames(corrL.PBMC[[jj]])){
      CorrTable.GeneX[jj, match(colnames(corrL.PBMC[[jj]]), genes.filtered)] = corrL.PBMC[[jj]] [, match(genes.filtered[ii], colnames(corrL.PBMC[[jj]]))]
    }
  }
  
  # Only include correlations found in more than 4 datasets
  index = which(colSums(!is.na(CorrTable.GeneX)) < 4)
  CorrTable.GeneX[,index] = NA
  
  #Calculate statistics for the gene 
  mean.x = colMeans(CorrTable.GeneX, na.rm=T)
  n.x = colSums(!is.na(CorrTable.GeneX),na.rm=T)
  sd.x = apply(FUN=sd, X=CorrTable.GeneX, MARGIN=2, na.rm=T)
  t.vector = mean.x*sqrt(n.x)/sd.x
  
  #Save t-values in matrix
  corrMx.PBMC[ii,] <- t.vector
  
}

# Save matrix to .csv file
setwd(myBigDataPath)
write.big.matrix(corrMx.PBMC, "corrMx_PBMC.csv", row.names = TRUE, col.names = TRUE, sep=',')


#### --------------------- Statistics -------------------------------------------
## Plot all gene-gene correlations in blood and PBMC
blood = as.matrix(corrMx.blood)
PBMC = as.matrix(corrMx.PBMC)

plot(blood[lower.tri(blood)], PBMC[lower.tri(PBMC)], main = "T-value for all gene-gene correlations", xlab = "t-value in blood", ylab = "t-value in PBMC")

## ---- Find interactions that only exist in one of the tissues  ----------------

# Only in PBMC -----------------------------------------------------------
CorrM1 = matrix(nrow = length(genes.filtered), ncol = length(genes.filtered))
rownames(CorrM1) = genes.filtered
colnames(CorrM1) = genes.filtered
for(i in 1:length(genes.filtered)){
  CorrM1[i,] = corrMx.PBMC[i,]-corrMx.blood[i,]
}

N = 5 #first N interactions
CorrL1 = list()
for(i in 1:ncol(CorrM1)){
  vx = CorrM1[-i,i]
  vx2 = sort(vx, decreasing=T)
  CorrL1[[i]] = rbind(rep(rownames(CorrM1)[i],N),names(vx2[1:N]),vx2[1:N])
}
TopCorrM1 = matrix(unlist(CorrL1),ncol=3,byrow=T)
colnames(TopCorrM1) = c("Gene 1", "Gene 2", "T-value")

# Plot distribution of t-values 
hist(as.numeric(TopCorrM1[,3]), breaks = length(TopCorrM1[,3]))


# Only in blood ------------------------------------------------------------
CorrM2 = matrix(nrow = length(genes.filtered), ncol = length(genes.filtered))
rownames(CorrM2) = genes.filtered
colnames(CorrM2) = genes.filtered
for(i in 1:length(genes.filtered)){
  CorrM2[i,] = corrMx.blood[i,]-corrMx.PBMC[i,]
}

N = 5 #first N interactions
CorrL2 = list()
for(i in 1:ncol(CorrM2)){
  vx = CorrM2[-i,i]
  vx2 = sort(vx, decreasing=T)
  CorrL2[[i]] = rbind(rep(rownames(CorrM2)[i],N),names(vx2[1:N]),vx2[1:N])
}
TopCorrM2 = matrix(unlist(CorrL2),ncol=3,byrow=T)
colnames(TopCorrM2) = c("Gene 1", "Gene 2", "T-value")

# Plot distribution of t-values
hist(as.numeric(TopCorrM2[,3]), breaks = length(TopCorrM2[,3]))


## Filter ToppCorrM  -------------------------------------

cutoff = 1 # Also tried with 0,6
remove = c()
for(j in 1:length(TopCorrM1[,3])){
  if(as.numeric(TopCorrM1[j,3]) < cutoff || is.na(TopCorrM1[j,3]) ){
    remove = append(remove, j)
  } else{
    next
  }
}
TopCorrM1_filtered = TopCorrM1[-remove,]

# Save table to file
write.csv(TopCorrM1_filtered, file="Supl.Table2_TopCorrPBMC.csv", sep = ',', col.names = TRUE)
#write.csv(TopCorrM1_filtered, file="Supl.Table4_TopCorrPBMC_cutoff0.6.csv", sep = ',', col.names = TRUE)


cutoff = 1 # Also tried with 0,6
remove = c()
for(j in 1:length(TopCorrM2[,3])){
  if(as.numeric(TopCorrM2[j,3]) < cutoff || is.na(TopCorrM2[j,3])){
    remove = append(remove, j)
  } else{
    next
  }
}
TopCorrM2_filtered = TopCorrM2[-remove,]

# Save table to file
write.csv(TopCorrM2_filtered, file="Supl.Table3_TopCorrBlood.csv", sep = ',', col.names = TRUE)
#write.csv(TopCorrM2_filtered, file="Supl.Table5_TopCorrBlood_cutoff0.6.csv", sep = ',', col.names = TRUE)


## -------------------- Create networks -------------------------------------

Corr.Graph1 = simplify(graph.data.frame(TopCorrM1_filtered, directed=F), edge.attr.comb="max")
Corr.Graph2 = simplify(graph.data.frame(TopCorrM2_filtered, directed=F), edge.attr.comb="max")

# Plot settings
V(Corr.Graph1)$label = ann.v[V(Corr.Graph1)$name]
V(Corr.Graph1)$label.family="sans" # type of font
V(Corr.Graph1)$label.cex=0.4 # size of font
V(Corr.Graph1)$label.dist = 0.5
V(Corr.Graph1)$label.degree = -pi/2
plot.igraph(Corr.Graph1, vertex.size=3, edge.arrow.size = 10, main = "Correlations in PBMC", sub = "N=5, cutoff=1" )


# Plot settings
V(Corr.Graph2)$label = ann.v[V(Corr.Graph2)$name]
V(Corr.Graph2)$label.family="sans" # type of font
V(Corr.Graph2)$label.cex=0.4 # size of font
V(Corr.Graph2)$label.dist = 0.5
V(Corr.Graph2)$label.degree = -pi/2
plot.igraph(Corr.Graph2, vertex.size=3, edge.arrow.size = 10, main = "Correlations in blood", sub = "N=5, cutoff=1" )


## --------------------------- Analysis ----------------------------------------------------------------

# What kind of genes are hub-genes? i.e genes with many interactions
top_genes1 = sort(degree(Corr.Graph1),decreasing=T)[1:20]
top_genes2 = sort(degree(Corr.Graph2),decreasing=T)[1:20]

## Summarize correlations for top 5 hub genes ------------------------------------------
top_genes1_table = matrix(nrow = 0, ncol = 2)
colnames(top_genes1_table) = c("Top gene", "Correlation gene")
for(i in 1:5) {
  for(j in 1:nrow(TopCorrM1_filtered)) {
    if (names(top_genes1[i]) == TopCorrM1_filtered[j, 1]) {
      top_genes1_table = rbind(top_genes1_table, c(names(top_genes1[i]),TopCorrM1_filtered[j, 2])) 
    } else if (names(top_genes1[i]) == TopCorrM1_filtered[j, 2]) {
      top_genes1_table = rbind(top_genes1_table, c(names(top_genes1[i]),TopCorrM1_filtered[j, 1])) 
    } else {
      next
    }
  }
}
top_genes1_table = top_genes1_table[!duplicated(top_genes1_table),]

# Save correlations for hub genes in .csv file 
write.csv(top_genes1_table, file="Supl.Table6_CorrelationsTop5HubGenes_PBMC.csv", sep = ',', col.names = TRUE)

top_genes2_table = matrix(nrow = 0, ncol = 2)
colnames(top_genes2_table) = c("Top gene", "Correlation gene")
for(i in 1:5) {
  for(j in 1:nrow(TopCorrM2_filtered)) {
    if (names(top_genes2[i]) == TopCorrM2_filtered[j, 1]) {
      top_genes2_table = rbind(top_genes2_table, c(names(top_genes2[i]),TopCorrM2_filtered[j, 2])) 
    } else if (names(top_genes2[i]) == TopCorrM2_filtered[j, 2]) {
      top_genes2_table = rbind(top_genes2_table, c(names(top_genes2[i]),TopCorrM2_filtered[j, 1])) 
    } else {
      next
    }
  }
}
top_genes2_table = top_genes2_table[!duplicated(top_genes2_table),]

# Save correlations for hub genes in .csv file 
write.csv(top_genes2_table, file="Supl.Table7_CorrelationsTop5HubGenes_Blood.csv", sep = ',', col.names = TRUE)


# ----------------- Extra ----------------

# In case R has crached:

# Read corrMx_blood and corrMx_PBMC
library(bigmemory)
setwd(myBigDataPath)
corrMx.blood = read.big.matrix("corrMx_blood.csv", sep = ",", header = TRUE, has.row.names = TRUE, descriptorfile = "corrMx_blood_descriptor.desc", backingfile = "corrMx_blood_backing.bin")
options(bigmemory.allow.dimnames=TRUE)
colnames(corrMx.blood) = genes.filtered
rownames(corrMx.blood) = genes.filtered
corrMx.PBMC = read.big.matrix("corrMx_PBMC.csv", sep = ",", header = FALSE, has.row.names = FALSE, descriptorfile = "corrMx_PBMC_descriptor.desc", backingfile = "corrMx_PBMC_backing.bin")
colnames(corrMx.PBMC) = genes.filtered
rownames(corrMx.PBMC) = genes.filtered


