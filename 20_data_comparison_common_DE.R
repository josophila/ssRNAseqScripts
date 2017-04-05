#This script is used to take lists of significant genes from multiple sets and detect
#whether any are in common.  In the current form of the data, up vs. down reg not taken
#into effect, which I want to do at some point -- 1/5/17

#Conversion to version 3.3 of genome. 
#Import list of genes with version 3.3 and corresponding version 1.6 gene names

rawsyn = read.csv('C:/Users/Josophila/Desktop/00_PhD/11_bioinformatics_tools/000_moss_genome/Ppatens/annotation/PhytozomeV11_download/annotation/syn_simple_prog.txt', 
                  header = FALSE, sep = '\t', stringsAsFactors = FALSE)
colnames(rawsyn) = c("version3_ID", "version1.6_ID")


#Import stress response data
coldStress = read.csv('C:/Users/Josophila/Desktop/SEQ_data/moss_omics_datasets/Cold_stress.txt', 
                      sep = '\t', header = TRUE, na.strings = '')

abaStress = read.csv('C:/Users/Josophila/Desktop/SEQ_data/moss_omics_datasets/aba_sheet.txt',
                     sep = '\t', header = TRUE, na.strings = '')

droughtStress = read.csv('C:/Users/Josophila/Desktop/SEQ_data/moss_omics_datasets/drought_stress.txt',
                         sep = '\t', header = TRUE, na.strings = '')

allstress = read.csv('C:/Users/Josophila/Desktop/SEQ_data/moss_omics_datasets/drought_stress.txt', 
                     sep = '\t', header = TRUE, na.strings = '')

#import gene atlas data
tis = c('spore', 'caul','chlor','gam','rhiz')
spore_sig = read.csv('./gene_atlas_data/geneAtlas_spore_significant_gen.txt', header = TRUE)
caul_sig = read.csv('./gene_atlas_data/geneAtlas_caul_significant_gen.txt', header = TRUE)
chlor_sig = read.csv('./gene_atlas_data/geneAtlas_chlor_significant_gen.txt', header = TRUE)
gam_sig = read.csv('./gene_atlas_data/geneAtlas_gam_significant_gen.txt', header = TRUE)
rhiz_sig = read.csv('./gene_atlas_data/geneAtlas_rhiz_significant_gen.txt', header = TRUE)


#makes a list, containing vectors that give the indices in rawsyn (which rows) have geneIDs 
#that are found in each column of the cold stress table.

sigGeneList = list(cold = coldStress, aba = abaStress, drought = droughtStress, 
                   all = allstress, spore = spore_sig, caulonema = caul_sig, 
                   chloronema = chlor_sig, gametophore = gam_sig, rhizoid = rhiz_sig)

#goes through the list of dataframes (sigGeneList), then compares contents of each column from dataframe to
#the synonym sheet.  This would actually be really easy to include more and more datasets on, since you can just
#tag them onto the thing that is now called 'sigGeneList'

namvec = c()        #empty vector to store names
maplist = list()    #empty list to store contents of match function (the map, essentially)
n = 0
for(j in 1:length(sigGeneList)){
  for(i in 1:ncol(sigGeneList[[j]])){
    n = n + 1
    namvec = c(namvec, paste(names(sigGeneList)[j], colnames(sigGeneList[[j]])[i], sep =''))
    maplist[[n]] = match(sigGeneList[[j]][,i], rawsyn$version1.6_ID)
  }
}
names(maplist) = namvec

#Do the same, detecting where in the synonym frame significant genes from each sc-rnaSeq
#are located. 

#step 1: import list of significant genes from each cluster into a list

clustdat = list()
clustnames = c()
for(i in 1:4){
  clustnames = c(clustnames, paste('cluster_', i, sep = ''))
  clustdat[[i]] = read.csv(paste('./sig_cluster_totals/c', i, '_6k_allsig.txt', sep = ''), header = FALSE)$V1
}
names(clustdat) = clustnames

#now map cluster reads from that list onto
cl_map = list()
for(i in 1:length(clustdat)){
  cl_map[[i]] = match(clustdat[[i]], rawsyn$version3_ID)
}
names(cl_map) = clustnames

totalNames = c(names(rawsyn), names(cl_map), names(maplist))
#PUT IT ALL TOGETHER BABY YEAH!
sig_compare = matrix(nrow = nrow(rawsyn), ncol = ncol(rawsyn) + length(cl_map) + length(maplist), NA)
sig_compare[,1] = rawsyn[,1]; sig_compare[,2] = rawsyn[,2]

colnames(sig_compare) = totalNames
for(i in 1:length(cl_map)){
  sig_compare[cl_map[[i]],i+2] = TRUE
}
for(i in 1:length(maplist)){
  sig_compare[maplist[[i]], i+6] = TRUE
}

#it's all together! let's write it, and open a new script to work on comparing cols
write.csv('./00_significant_by_set_raw.txt', row.names = FALSE, x = sig_compare, quote = FALSE)
