library(data.table)

JoesData <- read.csv("C:/Users/Josophila/Desktop/00_PhD/09_MHF_Transcriptomics/jan2014_annot_moss_reads-withV3.csv", header = TRUE, stringsAsFactors = FALSE)
JoesData <- data.frame(JoesData)

datasubset = subset(x = JoesData, select = c(1,3, 31:45)) # Subsets the data to include just gene names and RPM

##Next, want to take average of each tissue type, compare to average of all others (mimic analysis done by 10x)

#data for each group
leafdata = data.frame(datasubset$GeneID, datasubset[,3:5])
merdata = data.frame(datasubset$GeneID, datasubset$Ppa_Meristem_1.RPM, datasubset$Ppa_Meristem_2.RPM)
protdata = data.frame(datasubset$GeneID, datasubset$Ppa_Protonemal_Apical_Cell_1.RPM, datasubset$Ppa_Protonemal_Apical_Cell_2.RPM,datasubset$Ppa_Protonemal_Apical_Cell_3.RPM )
wholeplantdata = data.frame(datasubset$GeneID, datasubset$Ppa_Whole_Plant_1.RPM, datasubset$Ppa_Whole_Plant_2.RPM)
shootdata = data.frame (datasubset$GeneID, datasubset$Ppa_Shoot_1.RPM, datasubset$Ppa_Shoot_2.RPM, datasubset$Ppa_Shoot_3.RPM)
sporophytedata = data.frame(datasubset$GeneID, datasubset$Ppa_Sporophyte_2.RPM, datasubset$Ppa_Sporophyte_3.RPM)

#generate means for each group
leafmeans = data.frame(datasubset$GeneID, rowMeans(datasubset[,3:5]))
mermeans = data.frame(datasubset$GeneID, rowMeans(merdata[,2:3])) 
protmeans = data.frame(datasubset$GeneID, rowMeans(protdata[,2:4]))
wholeplantmeans = data.frame(datasubset$GeneID, rowMeans(wholeplantdata[,2:3]))
shootmeans = data.frame(datasubset$GeneID, rowMeans(shootdata[,2:4]))
spormeans = data.frame(datasubset$GeneID, rowMeans(sporophytedata[,2:3])) 

#data for compliments of each group
leafcompdata = data.frame(datasubset$GeneID, datasubset[,6:17])
mercompdata = data.frame(datasubset[,-c(6,7)])
protcompdata = data.frame(datasubset[,-c(8,9,10)])
wholeplantcompdata = data.frame(datasubset[,-c(11,12)])
shootcompdata = data.frame(datasubset[,-c(13,14,15)])
sporcompdata = data.frame(datasubset[,-c(16,17)])

#means of compliments of each group
leafcompmeans = data.frame(datasubset$GeneID, rowMeans(datasubset[,6:17]))
mercompmeans = data.frame(datasubset$GeneID, rowMeans(mercompdata[,-c(1,2)]))
protcompmeans = data.frame(datasubset$GeneID, rowMeans(protcompdata[,-c(1,2)]))
wholeplantcompmeans = data.frame(datasubset$GeneID, rowMeans(wholeplantcompdata[,-c(1,2)]))
shootcompmeans = data.frame(datasubset$GeneID, rowMeans(shootcompdata[,-c(1,2)]))
sporcompmeans = data.frame(datasubset$GeneID, rowMeans(sporcompdata[,-c(1,2)]))
#Test between samples
  #Make empty vectors to store pvalues
leafPData = c(); merPData = c(); protPData = c(); wholeplantPData = c(); shootPData = c(); sporPData = c()

#Loop through the pvalues for each set
for (i in seq(nrow(leafdata))){
  leafPData[i] = t.test(x = leafdata[i,2:4], y = leafcompdata[i,2:13])$p.value
  merPData[i] = t.test(x = merdata[i,-1], y = mercompdata[i, -c(1,2)])$p.value
  protPData[i] = t.test(x = protdata[i,-1], y = protcompdata[i, -c(1,2)])$p.value
  wholeplantPData[i] = t.test(x = wholeplantdata[i, -1], y = wholeplantcompdata[i, -c(1,2)])$p.value
  shootPData[i] = t.test(x = shootdata[i, -1], y = shootcompdata[i, -c(1,2)])$p.value
  sporPData[i] = t.test(x = sporophytedata[i, -1], y = sporcompdata[i, -c(1,2)])$p.value
}

#generate log2fc 
leafLog2 = c(); merlog2 = c(); protlog2 = c(); wholeplantlog2 = c(); shootlog2 = c(); sporlog2 = c()
for (i in seq(nrow(leafdata))){
  leafLog2[i] = log2(leafmeans[i,2]/leafcompmeans[i,2])
  merlog2[i] = log2(mermeans[i,2]/mercompmeans[i,2])
  protlog2[i] = log2(protmeans[i,2]/protcompmeans[i,2])
  wholeplantlog2[i] = log2(wholeplantmeans[i,2]/wholeplantcompmeans[i,2])
  shootlog2[i] = log2(shootmeans[i,2]/shootcompmeans[i,2])
  sporlog2[i] = log2(spormeans[i,2]/sporcompmeans[i,2])
}
#generate a dataframe with geneid, pvalue and log2fc 
names = JoesData$GeneID
leafallframe = data.frame(names, leafPData, leafLog2)
allframe = data.frame(names, leafPData, leafLog2, merPData, merlog2, protPData, protlog2, wholeplantPData, wholeplantlog2, shootPData, shootlog2, sporPData, sporlog2)

write.csv(allframe, file = "./mhf-10xstyle-comparison.csv")
