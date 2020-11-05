fn.x = "/Users/aina/Documents/X/Slutkurs/ITA/hsapiens.AllCodingHumanGenes.txt"
genes.all = sort(readLines(fn.x))

files = list.files(path = "/Users/aina/Documents/X/Slutkurs/ITA/Data/")

VarMatrix = matrix(nrow=length(files), ncol=length(genes.all))
colnames(VarMatrix) = genes.all 
rownames(VarMatrix) = files

setwd("/Users/aina/Documents/X/Slutkurs/ITA/Data/")
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
#Save variance matix to file
setwd("~/github/ITA/")
write.csv(VarMatrix, file ="VaranceMatrix.csv", row.names = TRUE)

# Find genes which needs to be removed from the matrix
bloodData = c("ExprM.GSE102459.2020_10_16.gz", "ExprM.GSE128224.2020_10_22.gz", "ExprM.GSE127792.2020_10_22.gz", "ExprM.GSE97590.2020_10_22.gz", "ExprM.GSE86627.2020_10_16.gz", "ExprM.GSE86331.2020_10_16.gz", "ExprM.GSE83951.2020_10_16.gz")
PBMCData = c("ExprM.GSE82152.2020_10_16.gz", "ExprM.GSE120115.2020_10_22.gz", "ExprM.GSE107990.2020_10_16.gz", "ExprM.GSE87186.2020_10_22.gz", "ExprM.GSE59743.2020_10_16.gz", "ExprM.GSE74816.2020_10_16.gz")

remove_2 = c()
for (i in 1:ncol(VarMatrix)){
  if(sum(is.na(VarMatrix[,i])) >=4){
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
    if( na_blood >=4 || na_PBMC >= 4){
      #remove = append(remove, colnames(VarMatrix)[i])
      remove_2 = append(remove_2, i)
    }
  }else{
    next
  }
}

#Remove genes from the matrix
modif_VarMatrix <- VarMatrix[,-remove_2]

#Calculate median for each column
VarMedian = c()
for(k in 1:ncol(modif_VarMatrix)){
  VarMedian = append(VarMedian, median(modif_VarMatrix[,k], na.rm=TRUE))
}
#View result in histogram 
hist(VarMedian, breaks = length(VarMedian), xlab = "Median with NA removed")



genes.filtered = genes.all[highvariance]

## Nu kan vi börja beräkna korrelationer på alla dataset och sammanställa resultaten.

cor.L = list()
for(i in 1:length(datasets.blood)){
	L1 = readLines(gzfile("ExprM.GSE107990.2020_10_16.gz")); closeAllConnections()
	L2 = strsplit(L1, split=",")
	Mx = sapply(FUN=function(v){return(as.numeric(v[-1]))},X=L2[2:length(L2)])
	rownames(Mx) = L2[[1]][-1]
	colnames(Mx) = sapply(FUN=getElement, X=L2[2:length(L2)], 1)
	for(i in 1:nrow(Mx)){Mx[i,] = Mx[i,]/max(Mx[i,],na.rm=T)}
	cor.L[[i]] = cor(Mx,use="pair")
}

## Sammanställ alla korrelationsvärden för en gen:

for(i in 1:length(genes.all)){

	
}