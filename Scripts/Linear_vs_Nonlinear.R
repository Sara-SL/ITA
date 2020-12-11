## Linear_vs_Nonlinear.R

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
genes.filtered = sort(readLines("genes.filtered_linear_vs_nonLinear.txt"))
PBMCData = c("ExprM.GSE82152.2020_10_16.gz", "ExprM.GSE120115.2020_10_22.gz", "ExprM.GSE87186.2020_10_22.gz", "ExprM.GSE59743.2020_10_16.gz", "ExprM.GSE74816.2020_10_16.gz")
datasets = PBMCData
geneNames = as.matrix(read.table("hsapiens.SYMBOL.txt", sep="\t", header=F))
ann.v = geneNames[,2]; 
names(ann.v) = geneNames[,1]

# For loading corrL.linear and corrL.nonLinear (if R has crashed) see secion "Extra"

#### -------------------- Create expression & correlation matrix for each dataset -------------------------------

## Linear correlation ---------------------------------------------------

# Create expression matrix
corrL.linear = list()
for(i in 1:length(datasets)){
  setwd(myDataPath)
  print(i)
  L1 = readLines(gzfile(datasets[i])); closeAllConnections()
  L2 = strsplit(L1, split=",")
  Mx = sapply(FUN=function(v){return(as.numeric(v[-1]))},X=L2[2:length(L2)])
  rownames(Mx) = L2[[1]][-1]
  colnames(Mx) = sapply(FUN=getElement, X=L2[2:length(L2)], 1)
  for(k in 1:nrow(Mx)){Mx[k,] = Mx[k,]/max(Mx[k,],na.rm=T)}
  
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
  corrL.linear = append(corrL.linear, list(cor(Mx, use="pairwise.complete.obs")))

  # Save correlation matrix to file in case R crashes 
  setwd(myBigDataPath)
  filename = paste0("corrL_linear_PBMC", i, ".csv")
  write.csv2(corrL.linear[i-1], filename)
}

# Remove correlation values below 0.7
for (i in 1:length(corrL.linear)) {
  print(i)
  corrL.linear[[i]] = replace(corrL.linear[[i]], corrL.linear[[i]] <= 0.7, NA)
}

# Initiate matrix for t-values
setwd(myBigDataPath)
corrMx.linear <- big.matrix(nrow = length(genes.filtered), ncol = length(genes.filtered), type = "double", init = NULL, backingfile = "corrMx_linear_backing.bin", descriptorfile = "corrMx_linear_descriptor.desc", shared = TRUE  )
options(bigmemory.allow.dimnames=TRUE)
colnames(corrMx.linear) = genes.filtered
rownames(corrMx.linear) = genes.filtered

# Create correlation table for each gene and calculate the t-value
for(ii in 1:length(genes.filtered)){
  print(ii)
  CorrTable.GeneX = matrix(nrow=length(datasets), ncol=length(genes.filtered))
  colnames(CorrTable.GeneX) = genes.filtered
  
  # Put genes in same order 
  for(jj in 1:length(corrL.linear)){ 
    if(genes.filtered[ii] %in% colnames(corrL.linear[[jj]])){
      CorrTable.GeneX[jj, match(colnames(corrL.linear[[jj]]), genes.filtered)] = corrL.linear[[jj]] [, match(genes.filtered[ii], colnames(corrL.linear[[jj]]))]
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
  corrMx.linear[ii,] <- t.vector
}

# Save matrix to .csv file
setwd(myBigDataPath)
write.big.matrix(corrMx.linear, "corrMx_linear.csv", row.names = TRUE, col.names = TRUE, sep=',')


## Non-linear correlation ------------------------------------

# Create expression matrix 
corrL.nonLinear = list()
for(i in 1:length(datasets)){
  setwd(myDataPath)
  print(i)
  L1 = readLines(gzfile(datasets[i])); closeAllConnections()
  L2 = strsplit(L1, split=",")
  Mx = sapply(FUN=function(v){return(as.numeric(v[-1]))},X=L2[2:length(L2)])
  rownames(Mx) = L2[[1]][-1]
  colnames(Mx) = sapply(FUN=getElement, X=L2[2:length(L2)], 1)
  for(k in 1:nrow(Mx)){Mx[k,] = Mx[k,]/max(Mx[k,],na.rm=T)}
  
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
  Mx_disc <- discretize(Mx, disc = "equalfreq")
  corrL.nonLinear.i = mutinformation(Mx_disc, method = "mm")
  
  # If all datsets could be run, then following line of code could save the resulting matrices in a list
  #corrL.nonLinear = append(corrL.nonLinear, list(corrL.nonLinear.i))
  # In our case the code in section "Extra" was used after running this forloop
  
  # Save correlation matrix to file in case R crashes 
  setwd(myBigDataPath)
  filename = paste0("corrL_nonLinear_PBMC", i, ".csv")
  write.csv2(corrL.nonLinear.i, filename)
}

# Initiate matrix for t-values 
setwd(myBigDataPath)
corrMx.nonLinear <- big.matrix(nrow = length(genes.filtered), ncol = length(genes.filtered), type = "double", init = NULL, backingfile = "corrMx_nonLinear_backing.bin", descriptorfile = "corrMx_nonLinear_descriptor.desc", shared = TRUE  )
options(bigmemory.allow.dimnames=TRUE)
colnames(corrMx.nonLinear) = genes.filtered
rownames(corrMx.nonLinear) = genes.filtered

# Create correlation table for each gene and calculate the t-value
for(ii in 1:length(genes.filtered)){
  print(ii)
  CorrTable.GeneX = matrix(nrow=length(datasets), ncol=length(genes.filtered))
  colnames(CorrTable.GeneX) = genes.filtered
  
  # Put genes in same order 
  for(jj in 1:length(corrL.nonLinear)){ 
    if(genes.filtered[ii] %in% colnames(corrL.nonLinear[[jj]])){
      CorrTable.GeneX[jj, match(colnames(corrL.nonLinear[[jj]]), genes.filtered)] = corrL.nonLinear[[jj]] [, match(genes.filtered[ii], colnames(corrL.nonLinear[[jj]]))]
    }
  }
  
  # Only include correlations found in more than 4 datasets
  index = which(colSums(!is.na(CorrTable.GeneX)) < 4)
  CorrTable.GeneX[,index] = NA
  
  #Calculate t-values for the gene 
  mean.x = colMeans(CorrTable.GeneX, na.rm=T)
  n.x = colSums(!is.na(CorrTable.GeneX),na.rm=T)
  sd.x = apply(FUN=sd, X=CorrTable.GeneX, MARGIN=2, na.rm=T)
  t.vector = mean.x*sqrt(n.x)/sd.x
  
  #Save t-values in matrix
  corrMx.nonLinear[ii,] <- t.vector
}

# Save matrix to .csv file
setwd(myBigDataPath)
write.big.matrix(corrMx.nonLinear, "corrMx_nonLinear.csv", row.names = TRUE, col.names = TRUE, sep=',')


#### --------------------- Statistics ----------------------------------
## Plot all gene-gene correlations in blood and PBMC
linear = as.matrix(corrMx.linear)
nonLinear = as.matrix(corrMx.nonLinear)

plot(linear[lower.tri(linear)], nonLinear[lower.tri(nonLinear)], ylim=c(0,40), main = "T-value for all gene-gene correlations", xlab = "t-value with Pearson", ylab = "t-value with MI")


# Only in linear -----------------------------------------------------------
CorrM1 = matrix(nrow = length(genes.filtered), ncol = length(genes.filtered))
rownames(CorrM1) = genes.filtered
colnames(CorrM1) = genes.filtered
for(i in 1:length(genes.filtered)){
  CorrM1[i,] = corrMx.linear[i,]-corrMx.nonLinear[i,]
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
hist(as.numeric(TopCorrM1[,3]), xlim = c(0,100), breaks = length(TopCorrM1[,3]), main = 'distribution of t-values with Pearson', xlab = 't-value')


# Only in non-linear ------------------------------------------------------------
CorrM2 = matrix(nrow = length(genes.filtered), ncol = length(genes.filtered))
rownames(CorrM2) = genes.filtered
colnames(CorrM2) = genes.filtered
for(i in 1:length(genes.filtered)){
  CorrM2[i,] = corrMx.nonLinear[i,]-corrMx.linear[i,]
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
hist(as.numeric(TopCorrM2[,3]), xlim = c(-150,50), breaks = length(TopCorrM2[,3]), main = 'distribution of t-values with MI', xlab = 't-value')


## Filter ToppCorrM  -------------------------------------

cutoff = 40
remove = c()
for(j in 1:length(TopCorrM1[,3])){
  if(as.numeric(TopCorrM1[j,3]) < cutoff || is.na(TopCorrM1[j,3]) ){
    remove = append(remove, j)
  } else{
    next
  }
}
TopCorrM1_filtered = TopCorrM1[-remove,]

cutoff = -40
remove = c()
for(j in 1:length(TopCorrM2[,3])){
  if(as.numeric(TopCorrM2[j,3]) < cutoff || is.na(TopCorrM2[j,3])){
    remove = append(remove, j)
  } else{
    next
  }
}
TopCorrM2_filtered = TopCorrM2[-remove,]


## -------------------- Create networks -------------------------------------## Create networks ------------------------------------------------------------

Corr.Graph1 = simplify(graph.data.frame(TopCorrM1_filtered, directed=F), edge.attr.comb="max")
Corr.Graph2 = simplify(graph.data.frame(TopCorrM2_filtered, directed=F), edge.attr.comb="max")


# Plot settings
V(Corr.Graph1)$label = ann.v[V(Corr.Graph1)$name]
V(Corr.Graph1)$label.family="sans" # type of font
V(Corr.Graph1)$label.cex=0.4 # size of font
V(Corr.Graph1)$label.dist = 0.5
V(Corr.Graph1)$label.degree = -pi/2
plot.igraph(Corr.Graph1, vertex.size=3, edge.arrow.size = 10, main = "Correlations with Pearsson", sub = "N=5, cutoff=40" )

# Plot settings
V(Corr.Graph2)$label = ann.v[V(Corr.Graph2)$name]
V(Corr.Graph2)$label.family="sans" # type of font
V(Corr.Graph2)$label.cex=0.4 # size of font
V(Corr.Graph2)$label.dist = 0.5
V(Corr.Graph2)$label.degree = -pi/2
plot.igraph(Corr.Graph2, vertex.size=3, edge.arrow.size = 10, main = "Correlations with MI", sub = "N=5, cutoff=-40" )


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
write.csv(top_genes1_table, file="CorrelationsTop5HubGenes_linear.csv", sep = ',', col.names = TRUE)


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
write.csv(top_genes2_table, file="CorrelationsTop5HubGenes_nonLinear.csv", sep = ',', col.names = TRUE)



# ----------------- Extra ----------------

# In case R has crached:

# Read corrL.linear files to a list
setwd(myBigDataPath)
corrL.linear = list()
corrL.linear = append(corrL.linear, list(read.csv2(file = 'corrL_linear_PBMC1.csv')))
corrL.linear = append(corrL.linear, list(read.csv2(file = 'corrL_linear_PBMC2.csv')))
corrL.linear = append(corrL.linear, list(read.csv2(file = 'corrL_linear_PBMC4.csv')))
corrL.linear = append(corrL.linear, list(read.csv2(file = 'corrL_linear_PBMC5.csv')))
corrL.linear = append(corrL.linear, list(read.csv2(file = 'corrL_linear_PBMC6.csv')))
for (i in 1:length(corrL.linear)) {
  rownames(corrL.linear[[i]]) = corrL.linear[[i]][,1]
  corrL.linear[[i]][1] <- NULL
}

# Read corrL.nonLinear files to a list
setwd(myBigDataPath)
corrL.nonLinear = list()
corrL.nonLinear = append(corrL.nonLinear, list(read.csv2(file = 'corrL_nonLinear_PBMC1.csv')))
corrL.nonLinear = append(corrL.nonLinear, list(read.csv2(file = 'corrL_nonLinear_PBMC2.csv')))
corrL.nonLinear = append(corrL.nonLinear, list(read.csv2(file = 'corrL_nonLinear_PBMC4.csv')))
corrL.nonLinear = append(corrL.nonLinear, list(read.csv2(file = 'corrL_nonLinear_PBMC5.csv')))
corrL.nonLinear = append(corrL.nonLinear, list(read.csv2(file = 'corrL_nonLinear_PBMC6.csv')))
for (i in 1:length(corrL.nonLinear)) {
  rownames(corrL.nonLinear[[i]]) = corrL.nonLinear[[i]][,1]
  corrL.nonLinear[[i]][1] <- NULL
}
