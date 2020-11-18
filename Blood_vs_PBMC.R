## Blood_vs_PBMC.R


## Set path ----------------------------------------------------------------

myPath = "~/github/ITA"
myDataPath = "~/github/ITA/Data/"
myBigDataPath = "~/github/ITA/BigData/"

## Load data ---------------------------------------------------------------

setwd(myPath)
genes.filtered = sort(readLines("genes.filtered.txt"))
bloodData = c("ExprM.GSE102459.2020_10_16.gz", "ExprM.GSE128224.2020_10_22.gz", "ExprM.GSE127792.2020_10_22.gz", "ExprM.GSE97590.2020_10_22.gz", "ExprM.GSE86627.2020_10_16.gz", "ExprM.GSE86331.2020_10_16.gz", "ExprM.GSE83951.2020_10_16.gz")
PBMCData = c("ExprM.GSE82152.2020_10_16.gz", "ExprM.GSE120115.2020_10_22.gz", "ExprM.GSE107990.2020_10_16.gz", "ExprM.GSE87186.2020_10_22.gz", "ExprM.GSE59743.2020_10_16.gz", "ExprM.GSE74816.2020_10_16.gz")
geneNames = as.matrix(read.table("hsapiens.SYMBOL.txt", sep="\t", header=F))
ann.v = geneNames[,2]; 
names(ann.v) = geneNames[,1]

# In case R crashes and we need to reload big matrices (saving environment doesn't seem to work)
install.packages("bigmemory")
library(bigmemory)
setwd(myBigDataPath)
corrMx.blood = read.big.matrix("corrMx_blood.csv", sep = ",", header = TRUE, has.row.names = TRUE, descriptorfile = "corrMx_blood_descriptor.desc", backingfile = "corrMx_blood_backing.bin")
options(bigmemory.allow.dimnames=TRUE)
colnames(corrMx.blood) = genes.filtered
rownames(corrMx.blood) = genes.filtered
corrMx.PBMC = read.big.matrix("corrMx_PBMC.csv", sep = ",", header = FALSE, has.row.names = FALSE, descriptorfile = "corrMx_PBMC_descriptor.desc", backingfile = "corrMx_PBMC_backing.bin")
colnames(corrMx.PBMC) = genes.filtered
rownames(corrMx.PBMC) = genes.filtered

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

# Remove correlation values below 0.7
for (i in 1:length(corrL.blood)) {
  print(i)
  corrL.blood[[i]] = replace(corrL.blood[[i]], corrL.blood[[i]] <= 0.7, NA)
}

# Initiate matrix for t-values
setwd(myBigDataPath)
corrMx.blood <- big.matrix(nrow = 14062, ncol = 14062, type = "double", init = NULL, dimnames = genes.filtered, backingfile = "corrMx_blood_backing.bin", descriptorfile = "corrMx_blood_descriptor.desc", shared = TRUE  )
options(bigmemory.allow.dimnames=TRUE)
colnames(corrMx.blood) = genes.filtered
rownames(corrMx.blood) = genes.filtered

# Create correlation table for each gene and calculate statistics
for(ii in 1:length(genes.filtered)){
  print(ii)
  CorrTable.GeneX = matrix(nrow=length(bloodData), ncol=length(genes.filtered))
  colnames(CorrTable.GeneX) = genes.filtered
  
  for(jj in 1:length(corrL.blood)){ 
    if(genes.filtered[ii] %in% colnames(corrL.blood[[jj]])){
      CorrTable.GeneX[jj, match(colnames(corrL.blood[[jj]]), genes.filtered)] = corrL.blood[[jj]] [, match(genes.filtered[ii], colnames(corrL.blood[[jj]]))]
    }
  }
  
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

# remove unessecary data to save memory 
rm(corrL.blood, L1, L2, Mx, remove_3, i, j, ii, jj, CorrTable.GeneX, mean.x, n.x, sd.x, t.vector)

# Save matrix to .csv file
setwd(myBigDataPath)
write.big.matrix(corrMx.blood, "corrMx_blood.csv", row.names = TRUE, col.names = TRUE, sep=',')


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

# Remove correlation values below 0.7
for (i in 1:length(corrL.PBMC)) {
  print(i)
  corrL.PBMC[[i]] = replace(corrL.PBMC[[i]], corrL.PBMC[[i]] <= 0.7, NA)
}

corrMx.PBMC <- big.matrix(nrow = 14062, ncol = 14062, type = "double", init = NULL, dimnames = genes.filtered, backingfile = "corrMx_PBMC_backing.bin", descriptorfile = "corrMx_PBMC_descriptor.desc", shared = TRUE )
options(bigmemory.allow.dimnames=TRUE)
colnames(corrMx.PBMC) = genes.filtered
rownames(corrMx.PBMC) = genes.filtered

# Create correlation table for each gene and calculate statistics
for(ii in 1:length(genes.filtered)){
  print(ii)
  CorrTable.GeneX = matrix(nrow=length(PBMCData), ncol=length(genes.filtered))
  colnames(CorrTable.GeneX) = genes.filtered
  
  for(jj in 1:length(corrL.PBMC)){ 
    if(genes.filtered[ii] %in% colnames(corrL.PBMC[[jj]])){
      CorrTable.GeneX[jj, match(colnames(corrL.PBMC[[jj]]), genes.filtered)] = corrL.PBMC[[jj]] [, match(genes.filtered[ii], colnames(corrL.PBMC[[jj]]))]
    }
  }
  
  index = which(colSums(!is.na(CorrTable.GeneX)) < 5)
  CorrTable.GeneX[,index] = NA
  
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

# Save matrix to .csv file
setwd(myBigDataPath)
write.big.matrix(corrMx.PBMC, "corrMx_PBMC.csv", row.names = TRUE, col.names = TRUE, sep=',')


# in 4 datasets -----------------------------------------------------------
# Initiate matrix for t-values
setwd(myBigDataPath)
corrMx.blood_2 <- big.matrix(nrow = 14062, ncol = 14062, type = "double", init = NULL, dimnames = genes.filtered, backingfile = "corrMx_blood_backing.bin", descriptorfile = "corrMx_blood_descriptor.desc", shared = TRUE  )
options(bigmemory.allow.dimnames=TRUE)
colnames(corrMx.blood_2) = genes.filtered
rownames(corrMx.blood_2) = genes.filtered

# Create correlation table for each gene and calculate statistics
for(ii in 1:length(genes.filtered)){
  print(ii)
  CorrTable.GeneX = matrix(nrow=length(bloodData), ncol=length(genes.filtered))
  colnames(CorrTable.GeneX) = genes.filtered
  
  for(jj in 1:length(corrL.blood)){ 
    if(genes.filtered[ii] %in% colnames(corrL.blood[[jj]])){
      CorrTable.GeneX[jj, match(colnames(corrL.blood[[jj]]), genes.filtered)] = corrL.blood[[jj]] [, match(genes.filtered[ii], colnames(corrL.blood[[jj]]))]
    }
  }
  
  index = which(colSums(!is.na(CorrTable.GeneX)) < 4)
  CorrTable.GeneX[,index] = NA
  
  #Calculate statistics for the gene 
  mean.x = colMeans(CorrTable.GeneX, na.rm=T)
  n.x = colSums(!is.na(CorrTable.GeneX),na.rm=T)
  sd.x = apply(FUN=sd, X=CorrTable.GeneX, MARGIN=2, na.rm=T)
  t.vector = mean.x*sqrt(n.x)/sd.x
  
  #Save t-values in matrix
  corrMx.blood_2[ii,] <- t.vector
}
--------------------------------------
corrMx.PBMC_2 <- big.matrix(nrow = 14062, ncol = 14062, type = "double", init = NULL, backingfile = "corrMx_PBMC2_backing.bin", descriptorfile = "corrMx_PBMC2_descriptor.desc", shared = TRUE )
options(bigmemory.allow.dimnames=TRUE)
colnames(corrMx.PBMC_2) = genes.filtered
rownames(corrMx.PBMC_2) = genes.filtered

# Create correlation table for each gene and calculate statistics
for(ii in 1:length(genes.filtered)){
  print(ii)
  CorrTable.GeneX = matrix(nrow=length(PBMCData), ncol=length(genes.filtered))
  colnames(CorrTable.GeneX) = genes.filtered
  
  for(jj in 1:length(corrL.PBMC)){ 
    if(genes.filtered[ii] %in% colnames(corrL.PBMC[[jj]])){
      CorrTable.GeneX[jj, match(colnames(corrL.PBMC[[jj]]), genes.filtered)] = corrL.PBMC[[jj]] [, match(genes.filtered[ii], colnames(corrL.PBMC[[jj]]))]
    }
  }
  
  index = which(colSums(!is.na(CorrTable.GeneX)) < 4)
  CorrTable.GeneX[,index] = NA
  
  #Calculate statistics for the gene 
  mean.x = colMeans(CorrTable.GeneX, na.rm=T)
  n.x = colSums(!is.na(CorrTable.GeneX),na.rm=T)
  sd.x = apply(FUN=sd, X=CorrTable.GeneX, MARGIN=2, na.rm=T)
  t.vector = mean.x*sqrt(n.x)/sd.x
  
  #Save t-values in matrix
  corrMx.PBMC_2[ii,] <- t.vector
  
}

### Statistics --------------------------------------------------------------

## Find interactions that only exist in one of the tissues  ----------------

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

#Distribution of t-values 
hist(as.numeric(TopCorrM1[,3]))
TopCorrM1_sort = TopCorrM1[order(TopCorrM1[,3], decreasing = T)]

# Only in blood ------------------------------------------------------------

CorrM2 = matrix(nrow = length(genes.filtered), ncol = length(genes.filtered))
rownames(CorrM2) = genes.filtered
colnames(CorrM2) = genes.filtered
for(i in 1:length(genes.filtered)){
  CorrM2[i,] = corrMx.blood[i,]-corrMx.PBMC[i,]
}

N = 5
CorrL2 = list()
for(i in 1:ncol(CorrM2)){
  vx = CorrM2[-i,i]
  vx2 = sort(vx, decreasing=T)
  CorrL2[[i]] = rbind(rep(rownames(CorrM2)[i],N),names(vx2[1:N]),vx2[1:N])
}
TopCorrM2 = matrix(unlist(CorrL2),ncol=3,byrow=T)
colnames(TopCorrM2) = c("Gene 1", "Gene 2", "T-value")

# Distribution of t-values
hist(as.numeric(TopCorrM2[,3]), xlim = c(0,100), breaks = length(TopCorrM2[,3]))


## Filter ToppCorrM based on histogram -------------------------------------

# Cutoff 
remove = c()
for(j in 1:length(TopCorrM1[,3])){
  if(as.numeric(TopCorrM1[j,3]) < 1 || is.na(TopCorrM1[j,3]) ){
    remove = append(remove, j)
  } else{
    next
  }
}
TopCorrM1_filtered = TopCorrM1[-remove,]

# Cutoff 
remove = c()
for(j in 1:length(TopCorrM2[,3])){
  if(as.numeric(TopCorrM2[j,3]) < 1 || is.na(TopCorrM2[j,3])){
    remove = append(remove, j)
  } else{
    next
  }
}
TopCorrM2_filtered = TopCorrM2[-remove,]

## Create graph ------------------------------------------------------------

install.packages("igraph")
library(igraph)

Corr.Graph1 = simplify(graph.data.frame(TopCorrM1_filtered, directed=F), edge.attr.comb="max")
Corr.Graph2 = simplify(graph.data.frame(TopCorrM2_filtered, directed=F), edge.attr.comb="max")


# Anta att gr.x är Ert nätverksobjekt från igraph.
V(Corr.Graph1)$label = ann.v[V(Corr.Graph1)$name]
V(Corr.Graph1)$label.family="sans" # Detta ger en lite snyggare font
V(Corr.Graph1)$label.cex=0.5 # Detta är storleken på fonten. Prova er fram tills ni hittar ett lämpligt värde.
V(Corr.Graph1)$label.dist = 1
plot.igraph(Corr.Graph1, vertex.size=3, main = "Correlations in PBMC", sub = "N=5, cutoff=1" )

# Anta att gr.x är Ert nätverksobjekt från igraph.
V(Corr.Graph2)$label = ann.v[V(Corr.Graph2)$name]
V(Corr.Graph2)$label.family="sans" # Detta ger en lite snyggare font
V(Corr.Graph2)$label.cex=0.5 # Detta är storleken på fonten. Prova er fram tills ni hittar ett lämpligt värde.
V(Corr.Graph2)$label.dist = 1
plot.igraph(Corr.Graph2, vertex.size=3, main = "Correlations in blood", sub = "N=5, cutoff=1" )

## Analysis ----------------------------------------------------------------

# What kind of genes are hub-genes? i.e genes with many interactions
top_genes1 = sort(degree(Corr.Graph1),decreasing=T)[1:20]
top_genes2 = sort(degree(Corr.Graph2),decreasing=T)[1:20]

# Measure of centrality based on shortest path
# how much a gene is affected by its surronding nodes
betweeness_genes1 = betweenness(Corr.Graph1, v = V(Corr.Graph1), directed = FALSE)
betweeness_genes2 = betweenness(Corr.Graph2, v = V(Corr.Graph2), directed = FALSE)

# Measures how many steps is required to access every other vertex from a given vertex.
# How independent a genes if from its surrounding genes
closenesse_genes1 = closeness(Corr.Graph1, vids = V(Corr.Graph1))
closenesse_genes2 = closeness(Corr.Graph2, vids = V(Corr.Graph2))

## Print results ------------------------------------------
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


# Extra  ------------------------------------------------------------------

commonCorr_table = matrix(nrow = 0, ncol = 2)
for (i in 1:nrow(TopCorrM1_filtered)) {
  for (j in 1: nrow(TopCorrM2_filtered)) {
    if (((TopCorrM1_filtered[i,1] == TopCorrM2_filtered[j,1]) && (TopCorrM1_filtered[i,2] == TopCorrM2_filtered[j,2]) ) || 
      ((TopCorrM1_filtered[i,2] == TopCorrM2_filtered[j,1]) && (TopCorrM1_filtered[i,1] == TopCorrM2_filtered[j,2]))) {
      commonCorr_table = rbind(commonCorr_table, TopCorrM1_filtered[i,1:2])
    }
  }
}


# Check max values of corrMx.blood and corrMx.PBMC ----------------------------------------------

max.blood = 0
for (a in 1:length(corrMx.blood[1,])) {
  print(a)
  if (is.na(max(corrMx.blood[,a])) ||max(corrMx.blood[,a]) == Inf ) {
    next
  }
  else if (max(corrMx.blood[,a]) >= max.blood ){
    max.blood = max(corrMx.blood[,a])
  } else {
    next
  }
}
# --> max.blood = 9904

max.PBMC = c()
for (a in 1:length(corrMx.PBMC[1,])) {
  print(a)
  corrMx.PBMC_2 = corrMx.PBMC[-a,a]
  if ( max(corrMx.PBMC_2[,a]) == Inf || is.na(max(corrMx.PBMC_2[,a]))) {
    next
  }
  else if (max(corrMx.PBMC_2[,a]) > 800 ){
    max.PBMC = append(max.PBMC, max(corrMx.PBMC_2[,a]))
    index = a
  } else {
    next
  }
}
# --> max.PBMC = 19134

CorrM3=corrMx.PBMC_2
N = 3 #first N interactions
CorrL3 = list()
for(i in 1:nrow(CorrM3)){
  vx = CorrM3[-i,i]
  vx2 = sort(vx, decreasing=T)
  CorrL3[[i]] = rbind(rep(rownames(CorrM3)[i],N),names(vx2[1:N]),vx2[1:N])
}
TopCorrM3 = matrix(unlist(CorrL3),ncol=3,byrow=T)
colnames(TopCorrM3) = c("Gene 1", "Gene 2", "T-value")

hist(as.numeric(TopCorrM3[,3]), xlim = c(0,100), ylim = c(0,100), breaks = length(TopCorrM3[,3]))

rm(vx, vx2, remove_3, N, min.PBMC, max.PBMC, L1, j, i,a, CorrL3, CorrM3, L2, Mx, TopCorrM3 )

for(i in 1:58) {
  corrMx.PBMC_2[i,i] = Inf
}

table = matrix(ncol = 3, nrow  = 6)
for (i in 1:length(PBMCData) ){
  table[i,1] = corrL.PBMC[[i]][match( "ENSG00000148965", rownames(corrL.PBMC[[i]]) ), match("ENSG00000235106", colnames(corrL.PBMC[[i]]) )]
  table[i,2] = corrL.PBMC[[i]][match( "ENSG00000162757", rownames(corrL.PBMC[[i]]) ), match("ENSG00000268043", colnames(corrL.PBMC[[i]]) )]
  table[i,3] = corrL.PBMC[[i]][match( "ENSG00000196369", rownames(corrL.PBMC[[i]]) ), match("ENSG00000204954", colnames(corrL.PBMC[[i]]) )]
}
