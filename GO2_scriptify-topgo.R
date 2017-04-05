#write Rscript to run through GO annotation for whatever input!
library("GO.db")
library("topGO")
library('Rgraphviz')

setwd('/Users/Josophila/Desktop/SEQ_data/10x_version3_genome_reanalysis/')

gotime = function(geneFile, outputFile, descript, set_type = "BP", numterms = 10){

sink(outputFile)

#Gene universe is physco version 3.3 from phytozome, from annotation file. 
geneID2GO = readMappings('./GO-map.txt', sep = '\t')
geneUniverse = names(geneID2GO)


#prepare data for inputting into TOPGO object.  
genesOfInterest = read.table(geneFile, header = FALSE)
genesOfInterest = as.character(genesOfInterest$V1)

geneList = factor(as.integer(geneUniverse %in% genesOfInterest))

names(geneList) = geneUniverse

#Build GO data object
# build the GOdata object in topGO
myGOdata <- new("topGOdata", description= descript, 
                ontology= set_type, allGenes=geneList,  
                annot = annFUN.gene2GO, gene2GO = geneID2GO)
myGOdata

# run the Fisher's exact tests
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")
resultElim <- runTest(myGOdata, algorithm="elim", statistic="fisher")
resultTopgo <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
resultParentchild <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")

# see how many results we get where weight01 gives a P-value <= 0.001:
mysummary <- summary(attributes(resultTopgo)$score <= 0.001)
numsignif <- as.integer(mysummary[[3]]) # how many terms is it true that P <= 0.001

# print out the top 'numsignif' results:
allRes <- GenTable(myGOdata, classicFisher = resultClassic, 
                   elimFisher = resultElim, topgoFisher = resultTopgo, 
                   parentchildFisher = resultParentchild, orderBy = "topgoFisher", 
                   ranksOf = "classicFisher", topNodes = numsignif)
allRes
topRes = GenTable(myGOdata, weightFisher = resultTopgo, orderBy = "resultTopGO", 
                   ranksOf = "weightFisher", topNodes = numterms)

# print a graph (to a pdf file) with the top 'numsignif' results:
output_file2 = paste(outputFile,"Topgo", sep="_")
printGraph(myGOdata, resultTopgo, firstSigNodes = numsignif, fn.prefix = output_file2, useInfo = "all", pdfSW = TRUE)

# print out the genes that are annotated with the significantly enriched GO terms:
myterms <- allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
for (i in 1:length(myterms))
{
  myterm <- myterms[i]
  mygenesforterm <- mygenes[myterm][[1]]
  myfactor <- mygenesforterm %in% genesOfInterest # find the genes that are in the list of genes of interest
  mygenesforterm2 <- mygenesforterm[myfactor == TRUE]
  mygenesforterm2 <- paste(mygenesforterm2, collapse=',')
  print(paste("Term",myterm,"genes:",mygenesforterm2))
}
# close the output file
sink() 

#write a file just giving the top howevermany, like 10, or 20, GO terms. The number is set by the 'numterms' variable.
write.csv(topRes, file = paste(outputFile, '.csv', sep=''), 
            row.names = FALSE)
}
