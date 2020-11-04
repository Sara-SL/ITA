fn.x = "~/github/ITA/hsapiens.AllCodingHumanGenes.txt"
genes.all = sort(readLines(fn.x))

files = list.files(path = "~/github/ITA/Data/")

VarMatrix = matrix(nrow=length(files), ncol=length(genes.all))
colnames(VarMatrix) = genes.all 
rownames(VarMatrix) = files

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

#
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