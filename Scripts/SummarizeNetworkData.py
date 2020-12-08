#### Script for summarizing data from networks

### Libraries to import
import csv
import pandas

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
file = open('/Users/saralundqvist/github/ITA/GeneData/hsapiens.SYMBOL.txt', 'r')
lines = file.readlines()
gene_names = []
for line in lines:
	split_line = line.strip().split("\t")
	gene_names.append(split_line)
file.close()

## Read files containing genes on PBMC with cutoff 1
file = open('/Users/saralundqvist/github/ITA/SupplementaryTables/Supl.Table2_TopCorrPBMC.csv', 'r')
csv_reader = csv.reader(file)
PBMC = list(csv_reader)
file.close()

## Read files containing genes on Blood with cutoff 1
file = open('/Users/saralundqvist/github/ITA/SupplementaryTables/Supl.Table3_TopCorrBlood.csv', 'r')
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
			PBMC_genes[row] = [PBMC_genes[row][0], gene_names[line][1], PBMC_genes[row][1]]

PBMC_genes.sort(key=lambda x:x[2], reverse=True)

# Write result to file
file = open("/Users/saralundqvist/github/ITA/gitignore/genes_in_PBMC.txt", "w")
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
			Blood_genes[row] = [Blood_genes[row][0], gene_names[line][1], Blood_genes[row][1]]

Blood_genes.sort(key=lambda x:x[2], reverse=True)

# Write result to file
file = open("/Users/saralundqvist/github/ITA/gitignore/genes_in_Blood.txt", "w")
csv.writer(file, delimiter=' ').writerows(Blood_genes)
file.close()


### Find common genes
common_genes = []
for i in range(len(PBMC_genes)):
	for j in range(len(Blood_genes)):
		if PBMC_genes[i][0] in Blood_genes[j]:
			common_genes.append([PBMC_genes[i][0], PBMC_genes[i][1]])

# Write to file
file = open("/Users/saralundqvist/github/ITA/gitignore/common_genes.txt", "w")
csv.writer(file, delimiter=' ').writerows(common_genes)
file.close()


### Genes specific for PBMC
specific_PBMC = []
for i in range(len(PBMC_genes)):
	if is_in(PBMC_genes[i][0], Blood_genes) == False:
		specific_PBMC.append([PBMC_genes[i][0], PBMC_genes[i][1]])

file = open("/Users/saralundqvist/github/ITA/gitignore/specific_PBMC.txt", "w")
csv.writer(file, delimiter=' ').writerows(specific_PBMC)
file.close()

### Genes specific for Blood
specific_Blood = []
for i in range(len(Blood_genes)):
	if is_in(Blood_genes[i][0], PBMC_genes) == False:
		specific_Blood.append([Blood_genes[i][0], Blood_genes[i][1]])

file = open("/Users/saralundqvist/github/ITA/gitignore/specific_Blood.txt", "w")
csv.writer(file, delimiter=' ').writerows(specific_Blood)
file.close()


# Create ToppCorrM1_filtered with gene names
topCorrM1_with_GeneNames = PBMC
for i in range(len(PBMC)):
	for j in range(len(gene_names)):
		if PBMC[i][1] == gene_names[j][0]:
			topCorrM1_with_GeneNames[i][1] = gene_names[j][1]
		elif PBMC[i][2] == gene_names[j][0]:
			topCorrM1_with_GeneNames[i][2] = gene_names[j][1]

file = open("/Users/saralundqvist/github/ITA/gitignore/TopCorrM1_filtered_with_GeneNames.csv", "w")
csv.writer(file, delimiter='\t').writerows(topCorrM1_with_GeneNames)
file.close()

# Create ToppCorrM2_filtered with gene names
topCorrM2_with_GeneNames = Blood
for i in range(len(Blood)):
	for j in range(len(gene_names)):
		if Blood[i][1] == gene_names[j][0]:
			topCorrM2_with_GeneNames[i][1] = gene_names[j][1]
		elif Blood[i][2] == gene_names[j][0]:
			topCorrM2_with_GeneNames[i][2] = gene_names[j][1]

file = open("/Users/saralundqvist/github/ITA/gitignore/TopCorrM2_filtered_with_GeneNames.csv", "w")
csv.writer(file, delimiter='\t').writerows(topCorrM2_with_GeneNames)
file.close()

#########  Create a summarizing table #########
#PBMC
topCorrPBMC = pandas.read_csv("/Users/saralundqvist/github/ITA/gitignore/TopCorrM1_filtered_with_GeneNames.csv", sep= '\t', header=None)
topCorrPBMC.columns = ['Tissue','Gene1', 'Gene2', 'T-value' ]
topCorrPBMC['Tissue'] = 'PBMC'
topCorrPBMC = topCorrPBMC.drop([0])

#Correlation count per gene
numbCorrPBMC = pandas.read_csv("/Users/saralundqvist/github/ITA/gitignore/genes_in_PBMC.txt", sep= ' ', header=None)
numbCorrPBMC.columns = ['ID','Gene', 'GeneCount']

topCorrPBMC['Gene1Count'] = ''
for i in range(1,len(topCorrPBMC['Gene1'])):
    for j in range(len(numbCorrPBMC['Gene'])):
        if (topCorrPBMC['Gene1'][i] == numbCorrPBMC['Gene'][j]):
            topCorrPBMC['Gene1Count'][i] = numbCorrPBMC['GeneCount'][j]

topCorrPBMC['Gene2Count'] = ''
for i in range(1,len(topCorrPBMC['Gene2'])):
    for j in range(len(numbCorrPBMC['Gene'])):
        if (topCorrPBMC['Gene2'][i] == numbCorrPBMC['Gene'][j]):
            topCorrPBMC['Gene2Count'][i] = numbCorrPBMC['GeneCount'][j]

#Common genes
commonGenes = pandas.read_csv("/Users/saralundqvist/github/ITA/gitignore/common_genes.txt", sep= ' ', header=None)
commonGenes.columns = ['ID','Gene']

topCorrPBMC['CommonGene1'] = ''
for i in range(1, len(topCorrPBMC['Gene1'])):
    for j in range(len(commonGenes['Gene'])):
        if topCorrPBMC['Gene1'][i] == commonGenes['Gene'][j]:
            topCorrPBMC['CommonGene1'][i] = 'YES'

topCorrPBMC['CommonGene2'] = ''
for i in range(1, len(topCorrPBMC['Gene2'])):
    for j in range(len(commonGenes['Gene'])):
        if topCorrPBMC['Gene2'][i] == commonGenes['Gene'][j]:
            topCorrPBMC['CommonGene2'][i] = 'YES'


#Blood
topCorrBlood = pandas.read_csv("/Users/saralundqvist/github/ITA/gitignore/TopCorrM2_filtered_with_GeneNames.csv", sep= '\t', header=None)
topCorrBlood.columns = ['Tissue','Gene1', 'Gene2', 'T-value' ]
topCorrBlood['Tissue'] = 'Blood'
topCorrBlood = topCorrBlood.drop([0])

#Correlation count per gene
numbCorrBlood = pandas.read_csv("/Users/saralundqvist/github/ITA/gitignore/genes_in_Blood.txt", sep= ' ', header=None)
numbCorrBlood.columns = ['ID','Gene', 'GeneCount']

topCorrBlood['Gene1Count'] = ''
for i in range(1,len(topCorrBlood['Gene1'])):
    for j in range(len(numbCorrBlood['Gene'])):
        if (topCorrBlood['Gene1'][i] == numbCorrBlood['Gene'][j]):
            topCorrBlood['Gene1Count'][i] = numbCorrBlood['GeneCount'][j]

topCorrBlood['Gene2Count'] = ''
for i in range(1,len(topCorrBlood['Gene2'])):
    for j in range(len(numbCorrBlood['Gene'])):
        if (topCorrBlood['Gene2'][i] == numbCorrBlood['Gene'][j]):
            topCorrBlood['Gene2Count'][i] = numbCorrBlood['GeneCount'][j]

#Common genes
commonGenes = pandas.read_csv("/Users/saralundqvist/github/ITA/gitignore/common_genes.txt", sep= ' ', header=None)
commonGenes.columns = ['ID','Gene']

topCorrBlood['CommonGene1'] = ''
for i in range(1, len(topCorrBlood['Gene1'])):
    for j in range(len(commonGenes['Gene'])):
        if topCorrBlood['Gene1'][i] == commonGenes['Gene'][j]:
            topCorrBlood['CommonGene1'][i] = 'YES'

topCorrBlood['CommonGene2'] = ''
for i in range(1, len(topCorrBlood['Gene2'])):
    for j in range(len(commonGenes['Gene'])):
        if topCorrBlood['Gene2'][i] == commonGenes['Gene'][j]:
            topCorrBlood['CommonGene2'][i] = 'YES'


#Merge DataFrames
sumTable = topCorrPBMC.append(topCorrBlood, ignore_index = True)
sumTable.to_csv(r'/Users/saralundqvist/github/ITA/SupplementaryTables/Supl.Table8_summarizingTable_blood_vs_PBMC.csv', sep='\t', index = False)



#########    .txt files to .csv files   #########
genesInPBMC = pandas.read_csv('/Users/saralundqvist/github/ITA/gitignore/genes_in_PBMC.txt', sep=" ", header=None)
genesInPBMC.columns = ['ID', 'Gene Name', 'Number of correlations']
genesInPBMC.to_csv(r'/Users/saralundqvist/github/ITA/SupplementaryTables/Supl.Table9_genes_in_PBMC.csv', sep='\t', index = False)

genesInBlood = pandas.read_csv('/Users/saralundqvist/github/ITA/gitignore/genes_in_Blood.txt', sep=" ", header=None)
genesInBlood.columns = ['ID', 'Gene Name', 'Number of correlations']
genesInBlood.to_csv(r'/Users/saralundqvist/github/ITA/SupplementaryTables/Supl.Table10_genes_in_Blood.csv', sep='\t', index = False)

commonGenes.columns = ['ID', 'Gene Name']
commonGenes.to_csv(r'/Users/saralundqvist/github/ITA/SupplementaryTables/Supl.Table11_common_genes.csv', sep='\t', index = False)

specificInPBMC = pandas.read_csv('/Users/saralundqvist/github/ITA/gitignore/specific_PBMC.txt', sep=" ", header=None)
specificInPBMC.columns = ['ID', 'Gene Name']
specificInPBMC.to_csv(r'/Users/saralundqvist/github/ITA/SupplementaryTables/Supl.Table12_specific_PBMC.csv', sep='\t', index = False)

specificInBlood = pandas.read_csv('/Users/saralundqvist/github/ITA/gitignore/specific_Blood.txt', sep=" ", header=None)
specificInBlood.columns = ['ID', 'Gene Name']
specificInBlood.to_csv(r'/Users/saralundqvist/github/ITA/SupplementaryTables/Supl.Table13_specific_Blood.csv', sep='\t', index = False)
