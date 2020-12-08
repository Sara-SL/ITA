## CorrelationTissue.R

#### ------------------ Set path ---------------------------
myPath = "/Users/aina/Documents/X/Slutkurs/github/ITA"
myDataPath = "/Users/aina/Documents/X/Slutkurs/github/ITA/Data/"
myBigDataPath = "/Users/aina/Documents/X/Slutkurs/github/ITA/myBigData/"
myResults = "/Users/aina/Documents/X/Slutkurs/github/ITA/Results"

#### ------------------ Load data -----------------------------
setwd(myPath)
genes.filtered = sort(readLines("/Users/aina/Documents/X/Slutkurs/github/ITA/genes.filtered_new.txt"))
bloodData = c("ExprM.GSE102459.2020_10_16.gz", "ExprM.GSE128224.2020_10_22.gz", "ExprM.GSE127792.2020_10_22.gz", "ExprM.GSE97590.2020_10_22.gz", "ExprM.GSE86627.2020_10_16.gz", "ExprM.GSE86331.2020_10_16.gz", "ExprM.GSE83951.2020_10_16.gz")
PBMCData = c("ExprM.GSE82152.2020_10_16.gz", "ExprM.GSE120115.2020_10_22.gz", "ExprM.GSE107990.2020_10_16.gz", "ExprM.GSE87186.2020_10_22.gz", "ExprM.GSE59743.2020_10_16.gz", "ExprM.GSE74816.2020_10_16.gz")
geneNames = as.matrix(read.table("hsapiens.SYMBOL.txt", sep="\t", header=F))
ann.v = geneNames[,2]; 
names(ann.v) = geneNames[,1]

#### In case R crashes and we need to reload big matrices (saving environment doesn't seem to work) -------
install.packages("bigmemory")
library(bigmemory)
setwd("/Users/aina/Documents/X/Slutkurs/github/ITA/Data/corrMx_without_NA")
corrMx.blood = read.big.matrix("corrMx_blood.csv", sep = ",", header = TRUE, has.row.names = TRUE, descriptorfile = "corrMx_blood_descriptor.desc", backingfile = "corrMx_blood_backing.bin")
options(bigmemory.allow.dimnames=TRUE)
colnames(corrMx.blood) = genes.filtered
rownames(corrMx.blood) = genes.filtered
corrMx.PBMC = read.big.matrix("corrMx_PBMC.csv", sep = ",", header = TRUE, has.row.names = TRUE, descriptorfile = "corrMx_PBMC_descriptor.desc", backingfile = "corrMx_PBMC_backing.bin")
colnames(corrMx.PBMC) = genes.filtered
rownames(corrMx.PBMC) = genes.filtered

#### -------------------- Create expression & correlation matrix for each dataset -------------------------------
### Blood -----------------------------------------------------------
setwd(myDataPath)
corrL.blood = list()
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
  #Only keep the filtered genes
  Mx <- Mx[, -remove_3]
  
  # Make a correlation matrix from the expression matrix
  corrL.blood = append(corrL.blood, list(cor(Mx, use="pairwise.complete.obs")))
}

# Without NA -----
# Used when NA is not used in the network
for(l in 1:length(corrL.blood)){
  corrL.blood[[l]] = replace(corrL.blood[[l]], corrL.blood[[l]]<=0.7, NA)
  print(l)
}

# With NA --------
# Set low correlating values to 0
for(l in 1:length(corrL.blood)){
  corrL.blood[[l]] = replace(corrL.blood[[l]], corrL.blood[[l]]<=0.7, 0)
}

## Initiate matrix for t-values -----
setwd(myBigDataPath)
corrMx.blood <- big.matrix(nrow = length(genes.filtered), ncol = length(genes.filtered), type = "double", init = NULL, dimnames = genes.filtered, backingfile = "corrMx_blood_backing.bin", descriptorfile = "corrMx_blood_descriptor.desc", shared = TRUE  )
options(bigmemory.allow.dimnames=TRUE)
colnames(corrMx.blood) = genes.filtered
rownames(corrMx.blood) = genes.filtered

# Create correlation table for each gene and calculate statistics
for(ii in 1:length(genes.filtered)){
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
  t.vector = c()
  for(i in 1:length(sd.x)){
    if(!is.na(sd.x[i]) && sd.x[i]== 0){
      t.vector =  append(t.vector, 0)
    }else{
      t.vector = append(t.vector, mean.x[i]*sqrt(n.x[i])/sd.x[i])
    }
  }
  
  #Save t-values in matrix
  corrMx.blood[ii,] <- t.vector
  
  print(ii)
}

# Save matrix to .csv file
setwd(myBigDataPath)
write.big.matrix(corrMx.blood, "corrMx_blood.csv", row.names = TRUE, col.names = TRUE, sep=',')

### PBMC -----------------------------------------------------------------------------------------
setwd(myDataPath)
corrL.PBMC = list()
for(i in 1:length(PBMCData)){
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

# Without NA ----
for(l in 1:length(corrL.PBMC)){
  corrL.PBMC[[l]] = replace(corrL.PBMC[[l]], corrL.PBMC[[l]]<=0.7, NA)
  print(l)
}

# With NA --------
#Set low correlating values to 0
for(l in 1:length(corrL.PBMC)){
  corrL.PBMC[[l]] = replace(corrL.PBMC[[l]], corrL.PBMC[[l]]<=0.7, 0)
  print(l)  
}

## Initiate matrix for t-values ------
setwd(myBigDataPath)
corrMx.PBMC <- big.matrix(nrow = length(genes.filtered), ncol = length(genes.filtered), type = "double", init = NULL, dimnames = genes.filtered, backingfile = "corrMx_PBMC_backing.bin", descriptorfile = "corrMx_PBMC_descriptor.desc", shared = TRUE  )
options(bigmemory.allow.dimnames=TRUE)
colnames(corrMx.PBMC) = genes.filtered
rownames(corrMx.PBMC) = genes.filtered

# Create correlation table for each gene and calculate statistics
for(ii in 1:length(genes.filtered)){
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
  for(i in 1:length(sd.x)){
    t.vector = c()
    if(!is.na(sd.x[i]) && sd.x == 0){
      t.vector = append(t.vector, 0)
    }else{
      t.vector = append(t.vector, mean.x[i]*sqrt(n.x[i])/sd.x[i])
    }
  }
  
  #Save t-values in matrix
  corrMx.PBMC[ii,] <- t.vector
  
  print(ii)
}

# Save matrix to .csv file
setwd(myBigDataPath)
write.big.matrix(corrMx.PBMC, "corrMx_PBMC.csv", row.names = TRUE, col.names = TRUE, sep=',')


#### --------------------- Statistics -------------------------------------------
### Find interactions that only exist in one of the tissues  ----------------

## Only in PBMC -----------------------------------------------------------
# Disregarding NA -----
CorrM1 = matrix(nrow = length(genes.filtered), ncol = length(genes.filtered))
rownames(CorrM1) = genes.filtered
colnames(CorrM1) = genes.filtered

for(i in 1:length(genes.filtered)){
  CorrM1[i,] = corrMx.PBMC[i,]-corrMx.blood[i,]
}

# Including NA ----
install.packages("hablar")
library(hablar)

CorrM1_NA = matrix(nrow = length(genes.filtered), ncol = length(genes.filtered))
rownames(CorrM1_NA) = genes.filtered
colnames(CorrM1_NA) = genes.filtered

for(i in 1:length(genes.filtered)){
  CorrM1_NA[i,] = corrMx.PBMC[i,] %minus_% corrMx.blood[i,]
  print(i)
}

# Put NA for the cells in CorrM1_NA which are included in NA_blood
for(i in 1:nrow(NA_blood)){
  CorrM1_NA[NA_blood[i,1], NA_blood[i,2]] = NA
}

# Replac inf with NA
CorrM1_NA = replace(CorrM1_NA, CorrM1_NA == Inf, NA)
CorrM1_NA = replace(CorrM1_NA, CorrM1_NA == -Inf, NA)

# Top correlations in PBMC disregarding NA ----
N = 5 
CorrL1 = list()
for(i in 1:ncol(CorrM1)){
  vx = CorrM1[-i,i]
  vx2 = sort(vx, decreasing=T)
  CorrL1[[i]] = rbind(rep(rownames(CorrM1)[i],N),names(vx2[1:N]),vx2[1:N])
}
TopCorrM1 = matrix(unlist(CorrL1),ncol=3,byrow=T)
colnames(TopCorrM1) = c("Gene 1", "Gene 2", "T-value")

#Distribution of t-values 
hist(as.numeric(TopCorrM1[,3]), main="Only PBMC disregarding NA")
TopCorrM1_sort = TopCorrM1[order(TopCorrM1[,3], decreasing = T)]

# Top correlations in PBMC including NA ----
N = 5 
CorrL1_NA = list()
for(i in 1:ncol(CorrM1_NA)){
  vx = CorrM1_NA[-i,i]
  vx2 = sort(vx, decreasing=T)
  CorrL1_NA[[i]] = rbind(rep(rownames(CorrM1_NA)[i],N),names(vx2[1:N]),vx2[1:N])
}
TopCorrM1_NA = matrix(unlist(CorrL1_NA),ncol=3,byrow=T)
colnames(TopCorrM1_NA) = c("Gene 1", "Gene 2", "T-value")

#Distribution of t-values 
hist(as.numeric(TopCorrM1_NA[,3]), main="Only PBMC including NA")
TopCorrM1_NA_sort = TopCorrM1_NA[order(TopCorrM1_NA[,3], decreasing = T)]

# Filter ToppCorrM1 based on histogram (disregarding NA) -------------------------------------
#Cutoff = 0,6
remove = c()
for(j in 1:length(TopCorrM1[,3])){
  if(as.numeric(TopCorrM1[j,3]) < 0.6 || is.na(TopCorrM1[j,3]) ){
    remove = append(remove, j)
  } else{
    next
  }
}
TopCorrM1_filtered_0.6 = TopCorrM1[-remove,]

#Cutoff = 1
remove = c()
for(j in 1:length(TopCorrM1[,3])){
  if(as.numeric(TopCorrM1[j,3]) < 1 || is.na(TopCorrM1[j,3]) ){
    remove = append(remove, j)
  } else{
    next
  }
}
TopCorrM1_filtered_1 = TopCorrM1[-remove,]


# Filter ToppCorrM1_NA based on histogram (including NA) ----------------------------
#Cutoff = 0,6
remove = c()
for(j in 1:length(TopCorrM1_NA[,3])){
  if(as.numeric(TopCorrM1_NA[j,3]) < 0.6 || is.na(TopCorrM1_NA[j,3]) ){
    remove = append(remove, j)
  } else{
    next
  }
}
TopCorrM1_NA_filtered_0.6 = TopCorrM1_NA[-remove,]

#Cutoff = 1
remove = c()
for(j in 1:length(TopCorrM1_NA[,3])){
  if(as.numeric(TopCorrM1_NA[j,3]) < 1 || is.na(TopCorrM1_NA[j,3]) ){
    remove = append(remove, j)
  } else{
    next
  }
}
TopCorrM1_NA_filtered_1 = TopCorrM1_NA[-remove,]

## Only in blood ------------------------------------------------------------
# Disregarding NA -----
CorrM2 = matrix(nrow = length(genes.filtered), ncol = length(genes.filtered))
rownames(CorrM2) = genes.filtered
colnames(CorrM2) = genes.filtered

for(i in 1:length(genes.filtered)){
  CorrM2[i,] = corrMx.blood[i,]-corrMx.PBMC[i,]
}

# Including NA -----
CorrM2_NA = matrix(nrow = length(genes.filtered), ncol = length(genes.filtered))
rownames(CorrM2_NA) = genes.filtered
colnames(CorrM2_NA) = genes.filtered

lowCorr.PBMC = matrix(ncol = 2, nrow = 0)
colnames(lowCorr.PBMC) = c("Gene_1", "Gene_2")

for(i in 1:length(genes.filtered)){
  for(j in 1:length(genes.filtered)){
    if(!is.na(corrMx.blood[i,j]) && !is.na(corrMx.PBMC[i,j]) && corrMx.PBMC[i,j] == 0 && corrMx.blood[i,j] != 0){
      CorrM2_NA[i,j] = corrMx.blood[i,j]
      lowCorr.PBMC = rbind(lowCorr.PBMC, c(colnames(corrMx.blood)[j], rownames(corrMx.blood)[i]))
    } else {
      CorrM2_NA[i,j] = corrMx.blood[i,j]-corrMx.PBMC[i,j]
    }
  }
}

# Top correlations in blood disregarding NA ----
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
hist(as.numeric(TopCorrM2[,3]), main = "Only Blood disregarding NA")

# Top correlations in blood including NA ----
N = 5
CorrL2_NA = list()
for(i in 1:ncol(CorrM2_NA)){
  vx = CorrM2_NA[-i,i]
  vx2 = sort(vx, decreasing=T)
  CorrL2_NA[[i]] = rbind(rep(rownames(CorrM2_NA)[i],N),names(vx2[1:N]),vx2[1:N])
}
TopCorrM2_NA = matrix(unlist(CorrL2_NA),ncol=3,byrow=T)
colnames(TopCorrM2_NA) = c("Gene 1", "Gene 2", "T-value")

# Distribution of t-values
hist(as.numeric(TopCorrM2_NA[,3]), main = "Only Blood including NA")

# Filter ToppCorrM2 based on histogram (disregarding NA) -------------------------------------
# Cutoff= 0.6
remove = c()
for(j in 1:length(TopCorrM2[,3])){
  if(as.numeric(TopCorrM2[j,3]) < 0.6 || is.na(TopCorrM2[j,3]) ){
    remove = append(remove, j)
  } else{
    next
  }
}
TopCorrM2_filtered_0.6 = TopCorrM2[-remove,]

# Cutoff = 1 
remove = c()
for(j in 1:length(TopCorrM2[,3])){
  if(as.numeric(TopCorrM2[j,3]) < 1 || is.na(TopCorrM2[j,3])){
    remove = append(remove, j)
  } else{
    next
  }
}
TopCorrM2_filtered_1 = TopCorrM2[-remove,]


# Filter ToppCorrM2 based on histogram (including NA) -------------------------------------
# Cutoff= 0.6
remove = c()
for(j in 1:length(TopCorrM2_NA[,3])){
  if(as.numeric(TopCorrM2_NA[j,3]) < 0.6 || is.na(TopCorrM2_NA[j,3]) ){
    remove = append(remove, j)
  } else{
    next
  }
}
TopCorrM2_NA_filtered_0.6 = TopCorrM2_NA[-remove,]

# Cutoff = 1 
remove = c()
for(j in 1:length(TopCorrM2_NA[,3])){
  if(as.numeric(TopCorrM2_NA[j,3]) < 1 || is.na(TopCorrM2_NA[j,3])){
    remove = append(remove, j)
  } else{
    next
  }
}
TopCorrM2_NA_filtered_1 = TopCorrM2_NA[-remove,]

### -------------------- Create graph -------------------------------------
install.packages("igraph")
library(igraph)

## ------ Graphs disregarding NA ------
# Graphs for cutoff 0.6 ----
Corr.Graph1_0.6 = simplify(graph.data.frame(TopCorrM1_filtered_0.6, directed=F), edge.attr.comb="max")
Corr.Graph2_0.6 = simplify(graph.data.frame(TopCorrM2_filtered_0.6, directed=F), edge.attr.comb="max")

V(Corr.Graph1_0.6)$label = ann.v[V(Corr.Graph1_0.6)$name]
V(Corr.Graph1_0.6)$label.family="sans" # Detta ger en lite snyggare font
V(Corr.Graph1_0.6)$label.cex=0.4 # Detta är storleken på fonten. Prova er fram tills ni hittar ett lämpligt värde.
V(Corr.Graph1_0.6)$label.dist = 0.5
V(Corr.Graph1_0.6)$label.degree = -pi/2
plot.igraph(Corr.Graph1_0.6, vertex.size=3, edge.arrow.size = 10, margin=0, main = "Correlations in PBMC", sub = "N=5, cutoff=0.6" )

V(Corr.Graph2_0.6)$label = ann.v[V(Corr.Graph2_0.6)$name]
V(Corr.Graph2_0.6)$label.family="sans" # Detta ger en lite snyggare font
V(Corr.Graph2_0.6)$label.cex=0.4 # Detta är storleken på fonten. Prova er fram tills ni hittar ett lämpligt värde.
V(Corr.Graph2_0.6)$label.dist = 0.5
V(Corr.Graph2_0.6)$label.degree = -pi/2
plot.igraph(Corr.Graph2_0.6, vertex.size=3, edge.arrow.size = 10, main = "Correlations in blood", sub = "N=5, cutoff=0.6" )

# Graphs for cutoff 1 ----
Corr.Graph1_1 = simplify(graph.data.frame(TopCorrM1_filtered_1, directed=F), edge.attr.comb="max")
Corr.Graph2_1 = simplify(graph.data.frame(TopCorrM2_filtered_1, directed=F), edge.attr.comb="max")

V(Corr.Graph1_1)$label = ann.v[V(Corr.Graph1_1)$name]
V(Corr.Graph1_1)$label.family="sans" # Detta ger en lite snyggare font
V(Corr.Graph1_1)$label.cex=0.4 # Detta är storleken på fonten. Prova er fram tills ni hittar ett lämpligt värde.
V(Corr.Graph1_1)$label.dist = 0.5
V(Corr.Graph1_1)$label.degree = -pi/2
plot.igraph(Corr.Graph1_1, vertex.size=3, edge.arrow.size = 10, main = "Correlations in PBMC", sub = "N=5, cutoff=1" )

V(Corr.Graph2_1)$label = ann.v[V(Corr.Graph2_1)$name]
V(Corr.Graph2_1)$label.family="sans" # Detta ger en lite snyggare font
V(Corr.Graph2_1)$label.cex=0.4 # Detta är storleken på fonten. Prova er fram tills ni hittar ett lämpligt värde.
V(Corr.Graph2_1)$label.dist = 0.5
V(Corr.Graph2_1)$label.degree = -pi/2
plot.igraph(Corr.Graph2_1, vertex.size=3, edge.arrow.size = 10, main = "Correlations in blood", sub = "N=5, cutoff=1" )


## ------- Graphs including NA ----------------
# Graphs for cutoff 0.6 ----
Corr.Graph1_NA_0.6 = simplify(graph.data.frame(TopCorrM1_NA_filtered_0.6, directed=F), edge.attr.comb="max")
Corr.Graph2_NA_0.6 = simplify(graph.data.frame(TopCorrM2_NA_filtered_0.6, directed=F), edge.attr.comb="max")

V(Corr.Graph1_NA_0.6)$label = ann.v[V(Corr.Graph1_NA_0.6)$name]
V(Corr.Graph1_NA_0.6)$label.family="sans" # Detta ger en lite snyggare font
V(Corr.Graph1_NA_0.6)$label.cex=0.4 # Detta är storleken på fonten. Prova er fram tills ni hittar ett lämpligt värde.
V(Corr.Graph1_NA_0.6)$label.dist = 0.5
V(Corr.Graph1_NA_0.6)$label.degree = -pi/2
plot.igraph(Corr.Graph1_NA_0.6, vertex.size=3, edge.arrow.size = 10, margin=0, main = "Correlations in PBMC including NA", sub = "N=5, cutoff=0.6" )

V(Corr.Graph2_NA_0.6)$label = ann.v[V(Corr.Graph2_NA_0.6)$name]
V(Corr.Graph2_NA_0.6)$label.family="sans" # Detta ger en lite snyggare font
V(Corr.Graph2_NA_0.6)$label.cex=0.4 # Detta är storleken på fonten. Prova er fram tills ni hittar ett lämpligt värde.
V(Corr.Graph2_NA_0.6)$label.dist = 0.5
V(Corr.Graph2_NA_0.6)$label.degree = -pi/2
plot.igraph(Corr.Graph2_NA_0.6, vertex.size=3, edge.arrow.size = 10, main = "Correlations in blood including NA", sub = "N=5, cutoff=0.6" )

# Graphs for cutoff 1 ----
Corr.Graph1_NA_1 = simplify(graph.data.frame(TopCorrM1_NA_filtered_1, directed=F), edge.attr.comb="max")
Corr.Graph2_NA_1 = simplify(graph.data.frame(TopCorrM2_NA_filtered_1, directed=F), edge.attr.comb="max")

V(Corr.Graph1_NA_1)$label = ann.v[V(Corr.Graph1_NA_1)$name]
V(Corr.Graph1_NA_1)$label.family="sans" # Detta ger en lite snyggare font
V(Corr.Graph1_NA_1)$label.cex=0.4 # Detta är storleken på fonten. Prova er fram tills ni hittar ett lämpligt värde.
V(Corr.Graph1_NA_1)$label.dist = 0.5
V(Corr.Graph1_NA_1)$label.degree = -pi/2
plot.igraph(Corr.Graph1_NA_1, vertex.size=3, edge.arrow.size = 10, main = "Correlations in PBMC including NA", sub = "N=5, cutoff=1" )

V(Corr.Graph2_NA_1)$label = ann.v[V(Corr.Graph2_NA_1)$name]
V(Corr.Graph2_NA_1)$label.family="sans" # Detta ger en lite snyggare font
V(Corr.Graph2_NA_1)$label.cex=0.4 # Detta är storleken på fonten. Prova er fram tills ni hittar ett lämpligt värde.
V(Corr.Graph2_NA_1)$label.dist = 0.5
V(Corr.Graph2_NA_1)$label.degree = -pi/2
plot.igraph(Corr.Graph2_NA_1, vertex.size=3, edge.arrow.size = 10, main = "Correlations in blood including NA", sub = "N=5, cutoff=1" )


## Analysis ----------------------------------------------------------------
# Top genes disregarding NA ----
# What kind of genes are hub-genes with cuttoff 0.6
top_genes1_0.6 = sort(degree(Corr.Graph1_0.6),decreasing=T)[1:20]
top_genes2_0.6 = sort(degree(Corr.Graph2_0.6),decreasing=T)[1:20]

# What kind of genes are hub-genes with cuttoff 1
top_genes1_1 = sort(degree(Corr.Graph1_1),decreasing=T)[1:20]
top_genes2_1 = sort(degree(Corr.Graph2_1),decreasing=T)[1:20]

# Top genes including NA ----
# What kind of genes are hub-genes with cuttoff 0.6
top_genes1_NA_0.6 = sort(degree(Corr.Graph1_NA_0.6),decreasing=T)[1:20]
top_genes2_NA_0.6 = sort(degree(Corr.Graph2_NA_0.6),decreasing=T)[1:20]

# What kind of genes are hub-genes with cuttoff 1
top_genes1_NA_1 = sort(degree(Corr.Graph1_NA_1),decreasing=T)[1:20]
top_genes2_NA_1 = sort(degree(Corr.Graph2_NA_1),decreasing=T)[1:20]






## ------- Get top gene interactions (disregarding NA) -----------

# PBMC - Blood cutoff 0.6----
top_genes1_table_0.6 = matrix(nrow = 0, ncol = 2)
colnames(top_genes1_table_0.6) = c("Top gene", "Correlation gene")
for(i in 1:5) {
  for(j in 1:nrow(TopCorrM1_filtered_0.6)) {
    if (names(top_genes1_0.6[i]) == TopCorrM1_filtered_0.6[j, 1]) {
      top_genes1_table_0.6 = rbind(top_genes1_table_0.6, c(names(top_genes1_0.6[i]),TopCorrM1_filtered_0.6[j, 2])) 
    } else if (names(top_genes1_0.6[i]) == TopCorrM1_filtered_0.6[j, 2]) {
      top_genes1_table_0.6 = rbind(top_genes1_table_0.6, c(names(top_genes1_0.6[i]),TopCorrM1_filtered_0.6[j, 1])) 
    } else {
      next
    }
  }
}
top_genes1_table_0.6 = top_genes1_table_0.6[!duplicated(top_genes1_table_0.6),]

write.csv(top_genes1_table_0.6, file="TopGeneCorr1_0.6.csv", sep = ',', col.names = TRUE)

# PBMC - Blood cutoff 1 ----
top_genes1_table_1 = matrix(nrow = 0, ncol = 2)
colnames(top_genes1_table_1) = c("Top gene", "Correlation gene")
for(i in 1:5) {
  for(j in 1:nrow(TopCorrM1_filtered_1)) {
    if (names(top_genes1_1[i]) == TopCorrM1_filtered_1[j, 1]) {
      top_genes1_table_1 = rbind(top_genes1_table_1, c(names(top_genes1_1[i]),TopCorrM1_filtered_1[j, 2])) 
    } else if (names(top_genes1_1[i]) == TopCorrM1_filtered_1[j, 2]) {
      top_genes1_table_1 = rbind(top_genes1_table_1, c(names(top_genes1_1[i]),TopCorrM1_filtered_1[j, 1])) 
    } else {
      next
    }
  }
}
top_genes1_table_1 = top_genes1_table_1[!duplicated(top_genes1_table_1),]

write.csv(top_genes1_table_1, file="TopGeneCorr1_1.csv", sep = ',', col.names = TRUE)


# Blood - PBMC cutoff 0.6 ----
top_genes2_table_0.6 = matrix(nrow = 0, ncol = 2)
colnames(top_genes2_table_0.6) = c("Top gene", "Correlation gene")
for(i in 1:5) {
  for(j in 1:nrow(TopCorrM2_filtered_0.6)) {
    if (names(top_genes2_0.6[i]) == TopCorrM2_filtered_0.6[j, 1]) {
      top_genes2_table_0.6 = rbind(top_genes2_table_0.6, c(names(top_genes2_0.6[i]),TopCorrM2_filtered_0.6[j, 2])) 
    } else if (names(top_genes2_0.6[i]) == TopCorrM2_filtered_0.6[j, 2]) {
      top_genes2_table_0.6 = rbind(top_genes2_table_0.6, c(names(top_genes2_0.6[i]),TopCorrM2_filtered_0.6[j, 1])) 
    } else {
      next
    }
  }
}
top_genes2_table_0.6 = top_genes2_table_0.6[!duplicated(top_genes2_table_0.6),]

write.csv(top_genes2_table_0.6, file="TopGeneCorr2_0.6.csv", sep = ',', col.names = TRUE)

# Blood - PBMC cutoff 1 ----
top_genes2_table_1 = matrix(nrow = 0, ncol = 2)
colnames(top_genes2_table_1) = c("Top gene", "Correlation gene")
for(i in 1:5) {
  for(j in 1:nrow(TopCorrM2_filtered_1)) {
    if (names(top_genes2_1[i]) == TopCorrM2_filtered_1[j, 1]) {
      top_genes2_table_1 = rbind(top_genes2_table_1, c(names(top_genes2_1[i]),TopCorrM2_filtered_1[j, 2])) 
    } else if (names(top_genes2_1[i]) == TopCorrM2_filtered_1[j, 2]) {
      top_genes2_table_1 = rbind(top_genes2_table_1, c(names(top_genes2_1[i]),TopCorrM2_filtered_1[j, 1])) 
    } else {
      next
    }
  }
}
top_genes2_table_1 = top_genes2_table_1[!duplicated(top_genes2_table_1),]

write.csv(top_genes2_table_1, file="TopGeneCorr2_1.csv", sep = ',', col.names = TRUE)


# Check for common gene-gene interactions -----
common_genes_1 = matrix(ncol=2, nrow=0)
colnames(common_genes_1) = c("Gene_1", "Gene_2")

for(i in 1:nrow(top_genes2_table_1)){
  for(j in 1:nrow(top_genes1_table_1)){
    if((top_genes2_table_1[i,1] == top_genes1_table_1[j,1] && top_genes2_table_1[i,2] == top_genes1_table_1[j,2]) ||
       (top_genes2_table_1[i,1] == top_genes1_table_1[j,2] && top_genes2_table_1[i,2] == top_genes1_table_1[j,1])){
      common_genes_1 = rbind(common_genes_1, top_genes1_table_1[j,1:2])
    } else {
      next
    }
  }
}

##### ---------- OLD -------------

# Analysis ----
# What kind of genes are hub-genes? i.e genes with many interactions
top_genes = sort(degree(Corr.Graph),decreasing=TRUE)[1:20]

'%notin%' <- Negate('%in%')
TopGeneTable = matrix(ncol = 3, nrow = 0)
for(i in 1:5){
  t_1 = TopCorrM_filtered[which(TopCorrM_filtered[,1] %in% names(top_genes[i])),]
  t_2 = TopCorrM_filtered[which(TopCorrM_filtered[,2] %in% names(top_genes[i])),]
  if(nrow(t_1) > nrow(t_2)){
    TopGeneTable = rbind(TopGeneTable, t_2)
    TopGeneTable = rbind(TopGeneTable, t_1[which(t_1[,2] %notin% t_2[,1]),])
  } else {
    TopGeneTable = rbind(TopGeneTable, t_1)
    TopGeneTable = rbind(TopGeneTable, t_2[which(t_2[,1] %notin% t_1[,2]),])
  }
}
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



# Find max in corrMx_blood
max = c()
for(a in 1:ncol(corrMx_test)){
  col = corrMx_test[-a,a]
  max = append(max, max(col))
}

# Convert NA to Inf
corrMx_test = corrMx
for(aa in 1:ncol(corrMx_test)){
  corrMx_test[aa,aa] = Inf
}


