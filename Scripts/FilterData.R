## FormatData

# Set path ----------------------------------------------------------------
myPath = "~/github/ITA/GeneData"
myDataPath = "~/github/ITA/ExpressionData/"
myTablePath = "~/github/ITA/Tables/"

# Load data ---------------------------------------------------------------
setwd(myPath)
fn.x = "hsapiens.AllCodingHumanGenes.txt"
genes.all = sort(readLines(fn.x))
files = list.files(path = myDataPath)


# Create varance matrix ---------------------------------------------------
VarMatrix = matrix(nrow=length(files), ncol=length(genes.all))
colnames(VarMatrix) = genes.all 
rownames(VarMatrix) = files

setwd(myDataPath)
for (i in 1:length(files)) {
  L1 = readLines(gzfile(files[i])); closeAllConnections()
  L2 = strsplit(L1, split=",")
  Mx = sapply(FUN=function(v){return(as.numeric(v[-1]))},X=L2[2:length(L2)])
  rownames(Mx) = L2[[1]][-1]
  colnames(Mx) = sapply(FUN=getElement, X=L2[2:length(L2)], 1)
  for(j in 1:nrow(Mx)) {Mx[j,] = Mx[j,]/max(Mx[j,],na.rm=T)}
  var.x = c()
  for(j in 1:ncol(Mx)){var.x[j] = var(Mx[,j])}
  names(var.x) = colnames(Mx)
  VarMatrix[i,] = var.x[genes.all]
}

# Save variance matix to file ---------------------------------------------
setwd(myTablePath)
write.csv(VarMatrix, file ="VaranceMatrix.csv", row.names = TRUE)


# Remove genes that are missing in too many datasets ----------------------
bloodData = c("ExprM.GSE102459.2020_10_16.gz", "ExprM.GSE128224.2020_10_22.gz", "ExprM.GSE127792.2020_10_22.gz", "ExprM.GSE97590.2020_10_22.gz", "ExprM.GSE86627.2020_10_16.gz", "ExprM.GSE86331.2020_10_16.gz", "ExprM.GSE83951.2020_10_16.gz")
PBMCData = c("ExprM.GSE82152.2020_10_16.gz", "ExprM.GSE120115.2020_10_22.gz", "ExprM.GSE107990.2020_10_16.gz", "ExprM.GSE87186.2020_10_22.gz", "ExprM.GSE59743.2020_10_16.gz", "ExprM.GSE74816.2020_10_16.gz")

# For blood vs PBMC analysis ----
# Each gene needs to exist in at least 4 datasets for each tissuetype 
remove= c()
for (i in 1:ncol(VarMatrix)){
  if(sum(is.na(VarMatrix[,i])) >2){
    na_blood = 0
    na_PBMC = 0
    for(j in 1:nrow(VarMatrix)){
      if(rownames(VarMatrix)[j] %in% bloodData && is.na(VarMatrix[j,i])){
        na_blood = na_blood +1
      }else if(rownames(VarMatrix)[j] %in% PBMCData && is.na(VarMatrix[j,i])){
        na_PBMC = na_PBMC+1
      }else{
        next
      }
    }
    if( na_blood >3 || na_PBMC >2){
      remove = append(remove, i)
    }
  }else{
    next
  }
}

# For linear vs non-linear correlation analysis ----
remove= c()
for (i in 1:ncol(VarMatrix)){
  if(sum(is.na(VarMatrix[,i])) > 0){
    remove = append(remove, i)
  }
}


# Remove genes from the matrix and genes.all----
modif_VarMatrix <- VarMatrix[,-remove]
genes.filtered = genes.all[-remove]


# Calculate median for each column(gene) ----------------------------------
VarMedian = c()
for(k in 1:ncol(modif_VarMatrix)){
  VarMedian = append(VarMedian, median(modif_VarMatrix[,k], na.rm=TRUE))
}

# View median distribution in histogram  ----------------------------------
hist(VarMedian, ylim=c(0,70), xlim=c(0,0.006), breaks = length(VarMedian), xlab = "Median variance with NA removed", main = 'Distribution of median variance')

# Remove genes that are at the bottom 25 or 40 percentile  ----------------------

#For blood vd PBMC analysis
percentile = quantile(VarMedian)

# For linear vs non-linear correlation analysis
percentile = quantile(VarMedian, probs = c(0, 0.4, 0.8, 1))

remove_2 = c()
for(l in 1:length(VarMedian)){
  if(VarMedian[l]<percentile[2]){
    remove_2 = append(remove_2, l)
  } else {
    next
  }
}
genes.filtered = genes.filtered[-remove_2]


# Save filtered genes to file ---------------------------------------------
setwd(myPath)
write.table(genes.filtered, file = "genes.filtered_blood_vs_PBMC.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
#or
write.table(genes.filtered, file = "genes.filtered_linear_vs_nonLinear.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
