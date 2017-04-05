##attempting to run topGO

library("GO.db")
library("topGO")
library('Rgraphviz')

#set working dir to one with v3 annotation info
setwd('/Users/Josophila/Desktop/SEQ_data/10x_version3_genome_reanalysis/')

#import genome annotation and other raw data 
ome = read.csv('./Ppatens_318_v3.3.annotation_info.txt',stringsAsFactors = FALSE, sep = '\t', header = TRUE)

#import singlecell cluster data
g3k = read.csv('./G04-3k_differential_expression.csv', sep = ',',header = TRUE)
# 
# #make just pacID and GO term, import mappings for topGO
# idgo = data.frame(ome$locusName, ome$GO)
# colnames(idgo) = c('id', 'GO')
# 
# write.table(file = './GO-map.txt', idgo, row.names = FALSE,
#           quote = FALSE, sep = '\t')


geneID2GO = readMappings('./GO-map.txt', sep = '\t')
geneUniverse = names(geneID2GO)

#cluster 1 genes of interest
genesOfInterest = data.frame(g3k$Gene.ID, g3k$Cluster.1.Adjusted.p.value)
genesOfInterest = subset(genesOfInterest, genesOfInterest$g3k.Cluster.1.Adjusted.p.value< .01)
genesOfInterest = genesOfInterest$g3k.Gene.ID
genesOfInterest = as.character(genesOfInterest)


#colnames(c1pval) = c('geneid', 'col1-p.val')
#IDlist = as.character(c1pval$geneid)

#make name and true/false (1/0) for whether gene detected in c1 has GO
geneList = factor(as.integer(geneUniverse%in%genesOfInterest))
names(geneList) = geneUniverse

#try to make this into a topGO data object
myGOdata = new("topGOdata", description = 'gransden3k cluster1' ,ontology = "BP", 
               allGenes = geneList,
               annot = annFUN.gene2GO , gene2GO = geneID2GO)
#yes!!!!

resultFisher = runTest(myGOdata, algorithm = 'weight01', statistic = 'fisher')

res10 = GenTable(myGOdata, weightFisher = resultFisher, orderBy = "resultFisher", 
                 ranksOf = "weightFisher", topNodes = 10)

showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 5, useInfo = 'all')

printGraph(myGOdata, resultFisher, firstSigNodes = 5, 
           fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
