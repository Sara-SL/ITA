#### Script for comparing the networks for PBMC and Blood 

### Libraries to import 
import csv

### Functions used in script
# Count how many times a given gene occurs in a list
def count_occurence(gene, list):
	count = 0
	
	for i in range(1,len(list)):
		for j in range(1,3):
			if list[i][j] == gene:
				count = count +1
	result = [gene , count]
	return result

# Checks if element in list of lists
def is_in(elem, list):
	result = False
	for i in range(len(list)):
		if elem in list[i]:
			result = True
	return result
	
# Find genes specific for one list
def specific(list1, list2):
	result = []
	for i in range(len(list1)):
		init = False
		for j in range(len(list2)):
			if list1[i][0] in list2[j]:
				init = True
		if init == False:
			result.append(list1[i][0])
	return result

### Read files  
## Read file with gene names 
file = open('/Users/aina/Documents/X/Slutkurs/github/ITA/hsapiens.SYMBOL.txt', 'r')
lines = file.readlines()	
gene_names = []
for line in lines:
	split_line = line.strip().split("\t")
	gene_names.append(split_line)
file.close()

## Read files containing genes on PBMC with cutoff 1
file = open('TopCorrM1_filtered_1.csv', 'r')
csv_reader = csv.reader(file)
PBMC = list(csv_reader)
file.close()

## Read files containing genes on Blood with cutoff 1
file = open('TopCorrM2_filtered_1.csv', 'r')
csv_reader = csv.reader(file)
Blood = list(csv_reader)
file.close()

### Extract which genes are present in respective network and number of interactions 
## For PBMC
PBMC_genes = []
for row in range(1,len(PBMC)):
	for colum in range(1,3): 
		if is_in(PBMC[row][colum],PBMC_genes) == False :
			line = count_occurence(PBMC[row][colum], PBMC)
			PBMC_genes.append(line)
			
for row in range(len(PBMC_genes)):
	for line in range(len(gene_names)):
		if PBMC_genes[row][0] == gene_names[line][0]:
			PBMC_genes[row][0] = gene_names[line][1]

PBMC_genes.sort(key=lambda x:x[1], reverse=True)

# Write result to file 
file = open("genes_in_PBMC.txt", "w")
csv.writer(file, delimiter=' ').writerows(PBMC_genes)
file.close()

## For Blood
Blood_genes = []
for row in range(1,len(Blood)):
	for column in range(1,3): 
		if is_in(Blood[row][colum],Blood_genes) == False :
			line = count_occurence(Blood[row][colum], Blood)
			Blood_genes.append(line)
			
for row in range(len(Blood_genes)):
	for line in range(len(gene_names)):
		if Blood_genes[row][0] == gene_names[line][0]:
			Blood_genes[row][0] = gene_names[line][1]

Blood_genes.sort(key=lambda x:x[1], reverse=True)		
	
# Write result to file 
file = open("genes_in_Blood.txt", "w")
csv.writer(file, delimiter=' ').writerows(Blood_genes)
file.close()

### Find common genes
common_genes = []
for i in range(len(PBMC_genes)):
	for j in range(len(Blood_genes)):
		if PBMC_genes[i][0] in Blood_genes[j]:
			common_genes.append(PBMC_genes[i][0])
# Write to file
file = open("common_gens.txt", "w")
for elem in common_genes:
	file.write(elem+'\n')
file.close()

### Genes specific for PBMC 
specific_PBMC = specific(PBMC_genes, Blood_genes)
file = open("specific_PBMC.txt", "w")
for elem in specific_PBMC:
	file.write(elem+'\n')
file.close()

### Genes specific for Blood 
specific_Blood = specific(Blood_genes, PBMC_genes)
file = open("specific_Blood.txt", "w")
for elem in specific_Blood:
	file.write(elem+'\n')
file.close()
