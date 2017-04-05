#preparation of data from gene atlas paper

#Import gene atlas data
genat = read.csv('C:/Users/Josophila/Desktop/SEQ_data/moss_omics_datasets/ppatens_atlas.txt',
                 sep = '\t', header = TRUE, stringsAsFactors = FALSE)

#take out some columns to just leave the data
unwanted = c('Mean', 'Calls','Enriched','Preferentially', 'SE',
             'Arquegonia', 'Sporophyte', 'Protein', 'GO_Term',
             'GeneName','Description', 'Ann_', 'Homolog','Iterative')
unwanted_map = c()
for(i in 1:length(unwanted)){
  unwanted_map = c(unwanted_map, grep(pattern = unwanted[i], colnames(genat)))
}

genDat = genat[,-unwanted_map]

#work with the gene atlas data somehow? What do I want to do with it...
#They have marked enriched as being higher than every other sample after doing pairwise comparisons
#seems pretty legit...  Should I do that with Margaret's data instead of what I did for 
#my reanalysis?  What I did was to pool all others, see if it is higher than all others. 
#I think that that's what 10x did...



tissList = c('Spores', 'Caulonema', 'Chloronema', 'Gametophore', 'Rhizoids')
tissMap = c()
spore_p =c(); caul_p = c(); chlor_p = c(); gam_p = c(); rhiz_p = c()
nameSafe = c()
for (i in 1:length(tissList)){
  foo_list = grep(pattern = tissList[i], colnames(genDat))
  foo = genDat[,foo_list]
  bar = genDat[,-c(1, foo_list)]
  for(j in 1:nrow(genDat)){
    if(i ==1){
      nameSafe = c(nameSafe, genDat[j,1])
      }
    testP = t.test(x = foo[j,],y=bar[j,])$p.value*length(tissList)   #adjusts p.value
    if(tissList[i] == 'Spores'){
      spore_p = c(spore_p, testP)
    }
    else if(tissList[i] == 'Caulonema'){
      caul_p = c(caul_p, testP)
    }
    else if(tissList[i] == 'Chloronema'){
      chlor_p = c(chlor_p, testP)
    }
    else if(tissList[i] == 'Gametophore'){
      gam_p = c(gam_p, testP)
    }
    else if(tissList[i] == 'Rhizoids'){
      rhiz_p = c(rhiz_p, testP)
    }
  }
}
geneatPRes = data.frame(nameSafe, spore_p, caul_p, chlor_p, gam_p, rhiz_p)


#write this CSV so it's saved.
write.csv(geneatPRes, file = './gene_atlas_data_p-adjusted.txt', 
          row.names = FALSE)

#get names of genes that are significant for each tissue type

tis = c('spore', 'caul','chlor','gam','rhiz')
sigAtList = list()

sigTissAll = function(tissue, pval = 0.001, frame){
  foo = paste(tissue, '_p', sep = '')
  subset(frame, frame[,foo] < pval)[,1]
}

for(i in 1:length(tis)){
  fnam = paste('./gene_atlas_data/geneAtlas_', tis[i], '_significant_gen.txt', sep = '')
  bloppy = sigTissAll(tissue = tis[i], frame = geneatPRes)
  write.csv(x = bloppy, file = fnam, row.names = FALSE, quote = FALSE)
}
