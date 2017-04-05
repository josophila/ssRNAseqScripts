#preparation of data from gene atlas paper
setwd('/Users/Josophila/Desktop/SEQ_data/10x_version3_genome_reanalysis/')

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
spore_p =c(); caul_p = c(); chlor_p = c(); gam_p = c(); rhiz_p = c()
spore_2fc =c(); caul_2fc = c(); chlor_2fc = c(); gam_2fc = c(); rhiz_2fc = c()

nameSafe = c()
for (i in 1:length(tissList)){
  foo_list = grep(pattern = tissList[i], colnames(genDat))    #this will look for the correct name in the headers
  foo = genDat[,foo_list]                                     #selects columns with correct name
  bar = genDat[,-c(1, foo_list)]                              #selects all the other columns
  for(j in 1:nrow(genDat)){
    if(i ==1){
      nameSafe = c(nameSafe, genDat[j,1])                     #takes the name from the right row, just to be consistent. prolly unnecessary
    }
    log2fc = log2(rowMeans(foo[j,]) / rowMeans(bar[j,]))
    testP = t.test(x = foo[j,],y=bar[j,])$p.value*length(tissList)   #*length(tissList) adjusts p.value 
    if(tissList[i] == 'Spores'){                      # setting log2f greater than 1.5 takes positively enriched only
      spore_p = c(spore_p, testP)
      spore_2fc = c(spore_2fc, log2fc)
    }
    else if(tissList[i] == 'Caulonema'){
      caul_p = c(caul_p, testP)
      caul_2fc = c(caul_2fc, log2fc)
    }
    else if(tissList[i] == 'Chloronema'){
      chlor_p = c(chlor_p, testP)
      chlor_2fc = c(chlor_2fc, log2fc)
    }
    else if(tissList[i] == 'Gametophore'){
      gam_p = c(gam_p, testP)
      gam_2fc = c(gam_2fc, log2fc)
    }
    else if(tissList[i] == 'Rhizoids'){
      rhiz_p = c(rhiz_p, testP)
      rhiz_2fc = c(rhiz_2fc, log2fc)
    }
  }
}
geneatPRes = data.frame(nameSafe, spore_p, spore_2fc,
                        caul_p, caul_2fc,
                        chlor_p, chlor_2fc,
                        gam_p, gam_2fc,
                        rhiz_p, rhiz_2fc)


#write this CSV so it's saved.
write.csv(geneatPRes, file = './outputs/11_gene_atlas_data_p-adjusted.txt', 
          row.names = FALSE)

#get names of genes that are significant for each tissue type.  sigTissAll is the function

sigTissAll = function(tissue, pval = 0.001, frame, posa = FALSE, posthresh = 1.5, nega = FALSE, negthresh = -1.5){
  ### Function subsets a data frame to include only data for which p is less than pthresh. Including pos or neg lets you extract only up or downregulated genes, given a default log2fc threshold of 1.5 or -1.5 to grant some robustness.
  
  #First, if pos and neg are set to true, return an error
  if(posa == TRUE & nega == TRUE){
    
    print("You can't go up AND down, silly. Set only pos or nega to true to select only up or down-regulated genes, respectively.  Otherwise, leave both to be false to select all differentially expressed (p < pval) genes regardless of their log2fc (relative to average of other genes in set)")
    
  }
  else{
  foo = paste(tissue, '_p', sep = '')
  woo = paste(tissue, '_2fc', sep = '')
  if(posa == TRUE){
    subset(frame, frame[,foo] < pval & frame[,woo] > posthresh)[,1]
    #print("it's running positive")
  }
  else if(nega == TRUE){
    subset(frame, frame[,foo] < pval & frame[,woo] < negthresh)[,1]
    #print("it's running negative")
  }
  else if (nega == FALSE & posa == FALSE){
    subset(frame, frame[,foo] < pval)[,1]
    #print("it's running without a log2fc filter")
    }
  }
}

#list of tissue names to iterate over
tis = c('spore', 'caul','chlor','gam','rhiz')

#get all sig
for(i in 1:length(tis)){
  fnam = paste('./outputs/11_gene_atlas_data/geneAtlas_', tis[i], '_significant_gen.txt', sep = '')
  bloppy = sigTissAll(tissue = tis[i], frame = geneatPRes, posa = FALSE, nega = FALSE)
  write.csv(x = bloppy, file = fnam, row.names = FALSE, quote = FALSE)
}

#get positive sig (genes in each tissue that are above average)
for(i in 1:length(tis)){
  fnam = paste('./outputs/11_gene_atlas_data/geneAtlas_', tis[i], '_significant_gen_UP.txt', sep = '')
  bloppy = sigTissAll(tissue = tis[i], frame = geneatPRes, posa = TRUE, nega = FALSE)
  write.csv(x = bloppy, file = fnam, row.names = FALSE, quote = FALSE)
}

#get negative (genes in each tissue that are below average)
for(i in 1:length(tis)){
  fnam = paste('./outputs/11_gene_atlas_data/geneAtlas_', tis[i], '_significant_gen_DOWN.txt', sep = '')
  bloppy = sigTissAll(tissue = tis[i], frame = geneatPRes, nega = TRUE)
  write.csv(x = bloppy, file = fnam, row.names = FALSE, quote = FALSE)
}
