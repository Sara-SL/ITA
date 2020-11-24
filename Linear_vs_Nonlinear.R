## Linear_vs_Nonlinear.R

library(infotheo)
library(bigmemory)

## Set path ----------------------------------------------------------------

myPath = "~/github/ITA"
myDataPath = "~/github/ITA/Data/"
myBigDataPath = "~/github/ITA/BigData/"

## Load data ---------------------------------------------------------------

setwd(myPath)
genes.filtered = sort(readLines("genes.filtered.txt"))
datasets = list.files(path = myDataPath)
geneNames = as.matrix(read.table("hsapiens.SYMBOL.txt", sep="\t", header=F))
ann.v = geneNames[,2]; 
names(ann.v) = geneNames[,1]

# Linear correlation ------------------------------

# Create expression matrix
setwd(myDataPath)
corrL.linear = list()
for(i in 1:length(datasets)){
  print(i)
  L1 = readLines(gzfile(datasets[i])); closeAllConnections()
  L2 = strsplit(L1, split=",")
  Mx = sapply(FUN=function(v){return(as.numeric(v[-1]))},X=L2[2:length(L2)])
  rownames(Mx) = L2[[1]][-1]
  colnames(Mx) = sapply(FUN=getElement, X=L2[2:length(L2)], 1)
  for(i in 1:nrow(Mx)){Mx[i,] = Mx[i,]/max(Mx[i,],na.rm=T)}
  
  #Only keep the filtered genes
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
  corrL.linear = append(corrL.linear, list(cor(Mx, use="pairwise.complete.obs")))
}

# Remove correlation values below 0.7
for (i in 1:length(corrL.linear)) {
  print(i)
  corrL.linear[[i]] = replace(corrL.linear[[i]], corrL.linear[[i]] <= 0.7, NA)
}

# Initiate matrix for t-values
setwd(myBigDataPath)
corrMx.linear <- big.matrix(nrow = 14062, ncol = 14062, type = "double", init = NULL, backingfile = "corrMx_linear_backing.bin", descriptorfile = "corrMx_linear_descriptor.desc", shared = TRUE  )
options(bigmemory.allow.dimnames=TRUE)
colnames(corrMx.linear) = genes.filtered
rownames(corrMx.linear) = genes.filtered

# Create correlation table for each gene and calculate statistics
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

rm(corrL.linear)
# Non-linear correlation ------------------------------------

# Create expression matrix 
setwd(myDataPath)
corrL.nonLinear = list()
for(i in 1:length(datasets)){
  print(i)
  L1 = readLines(gzfile(datasets[i])); closeAllConnections()
  L2 = strsplit(L1, split=",")
  Mx = sapply(FUN=function(v){return(as.numeric(v[-1]))},X=L2[2:length(L2)])
  rownames(Mx) = L2[[1]][-1]
  colnames(Mx) = sapply(FUN=getElement, X=L2[2:length(L2)], 1)
  for(i in 1:nrow(Mx)){Mx[i,] = Mx[i,]/max(Mx[i,],na.rm=T)}
  
  # Filter genes
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
  Mx_disc <- discretize(Mx, disc = "equalfreq")
  corrL.nonLinear.i = mutinformation(Mx_disc, method = "mm")
  corrL.nonLinear = append(corrL.nonLinear, list(corrL.nonLinear.i))
  
  # Save correlation matrix to file in case R crashes 
  filename = paste0("corrL_nonLinear_", i, ".csv")
  write.csv2(corrL.nonLinear.i, filename)
}

# Initiate matrix for t-values
setwd(myBigDataPath)
corrMx.nonLinear <- big.matrix(nrow = 14062, ncol = 14062, type = "double", init = NULL, backingfile = "corrMx_nonLinear_backing.bin", descriptorfile = "corrMx_nonLinear_descriptor.desc", shared = TRUE  )
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


### Statistics --------------------------------------------------------------

## Find interactions that only exist in one of the tissues  ----------------

# Only in PBMC -----------------------------------------------------------

CorrM1 = matrix(nrow = length(genes.filtered), ncol = length(genes.filtered))
rownames(CorrM1) = genes.filtered
colnames(CorrM1) = genes.filtered
for(i in 1:length(genes.filtered)){
  CorrM1[i,] = corrMx2.PBMC[i,]-corrMx2.blood[i,]
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
  CorrM2[i,] = corrMx2.blood[i,]-corrMx2.PBMC[i,]
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

N = 10
time = c()
for (i in 1:8) {
  print(i)
  Mx_test = Mx[,1:N]
  time_start = Sys.time()
  Mx_disc = discretize(Mx_test, disc = "equalfreq")
  corrL.nonLinear.i = mutinformation(Mx_disc, method = "mm")
  time_end = Sys.time()
  N = N*2
  time = c(time, time_end-time_start)
}

x = c(10,20,40,80,160,320,640,1280)
plot(x, time, xlab = "number of columns", ylab = "time")
plot(time)
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

