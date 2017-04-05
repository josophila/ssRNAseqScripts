#This script is used to take lists of significant genes from multiple sets and detect
#whether any are in common.  In the current form of the data, up vs. down reg not taken
#into effect, which I want to do at some point -- 1/5/17

#set working directory
setwd('/Users/Josophila/Desktop/SEQ_data/10x_version3_genome_reanalysis/')

#Conversion to version 3.3 of genome. 
#Import list of genes with version 3.3 and corresponding version 1.6 gene names

rawsyn = read.csv('C:/Users/Josophila/Desktop/00_PhD/11_bioinformatics_tools/000_moss_genome/Ppatens/annotation/PhytozomeV11_download/annotation/syn_simple_prog.txt', header = FALSE, sep = '\t', stringsAsFactors = FALSE)
colnames(rawsyn) = c("version3_ID", "version1.6_ID")


#Import stress response data
stressList = c('ABA', 'Cold', 'Drought', 'Salt')
tpts = c('.0.5.h', '.4.0.h')
updog = c('UP', 'DOWN')

#make a bunch of empty things with names!

aba05_up =read.csv('./outputs/stress_data/ABA.0.5.hUP.txt')
aba4_up =read.csv('./outputs/stress_data/ABA.4.0.hUP.txt')
aba05_down =read.csv('./outputs/stress_data/ABA.0.5.hDOWN.txt')
aba4_down = read.csv('./outputs/stress_data/ABA.4.0.hDOWN.txt')
cold05_up =read.csv('./outputs/stress_data/Cold.0.5.hUP.txt')
cold4_up =read.csv('./outputs/stress_data/Cold.4.0.hUP.txt')
cold05_down =read.csv('./outputs/stress_data/Cold.0.5.hDOWN.txt')
cold4_down =read.csv('./outputs/stress_data/Cold.4.0.hDOWN.txt')
drought05_up =read.csv('./outputs/stress_data/Drought.0.5.hUP.txt')
drough4_up =read.csv('./outputs/stress_data/Drought.4.0.hUP.txt') 
drought05_down =read.csv('./outputs/stress_data/Drought.0.5.hDOWN.txt')
drought4_down =read.csv('./outputs/stress_data/Drought.4.0.hDOWN.txt ')
salt05_up =read.csv('./outputs/stress_data/Salt.0.5.hUP.txt')
salt4_up =read.csv('./outputs/stress_data/Salt.4.0.hUP.txt')
salt05_down =read.csv('./outputs/stress_data/Salt.0.5.hDOWN.txt')
salt4_down =read.csv('./outputs/stress_data/Salt.4.0.hDOWN.txt ')

stressvec = c(aba05_up, aba4_up, aba05_down, aba4_down, cold05_up, cold4_up, cold05_down, cold4_down, drought05_up, drough4_up, drought05_down, drought4_down, salt05_up, salt4_up, salt05_down, salt4_down)


#import gene atlas data
tis = c('spore', 'caul','chlor','gam','rhiz')
spore_sig = read.csv('./outputs/11_gene_atlas_data/geneAtlas_spore_significant_gen.txt', header = TRUE)
caul_sig = read.csv('./outputs/11_gene_atlas_data/geneAtlas_caul_significant_gen.txt', header = TRUE)
chlor_sig = read.csv('./outputs/11_gene_atlas_data/geneAtlas_chlor_significant_gen.txt', header = TRUE)
gam_sig = read.csv('./outputs/11_gene_atlas_data/geneAtlas_gam_significant_gen.txt', header = TRUE)
rhiz_sig = read.csv('./outputs/11_gene_atlas_data/geneAtlas_rhiz_significant_gen.txt', header = TRUE)

spore_up = read.csv('./outputs/11_gene_atlas_data/geneAtlas_spore_significant_gen_UP.txt', header = TRUE)
caul_up = read.csv('./outputs/11_gene_atlas_data/geneAtlas_caul_significant_gen_UP.txt', header = TRUE)
chlor_up = read.csv('./outputs/11_gene_atlas_data/geneAtlas_chlor_significant_gen_UP.txt', header = TRUE)
gam_up = read.csv('./outputs/11_gene_atlas_data/geneAtlas_gam_significant_gen_UP.txt', header = TRUE)
rhiz_up = read.csv('./outputs/11_gene_atlas_data/geneAtlas_rhiz_significant_gen_UP.txt', header = TRUE)

#makes a list, containing vectors that give the indices in rawsyn (which rows) have geneIDs 
#that are found in each column of the cold stress table.

sigGeneList = list(aba_0.5hr_up = aba05_up, aba_4.0hr_up = aba4_up, aba_0.5hr_down =  aba05_down, aba_4.0hr_up = aba4_down,
                   cold_0.5hr_up = cold05_up, cold_4.0hr_up = cold4_up, cold_0.5hr_down = cold05_down, cold_4.0hr_down = cold4_down,
                   drought_0.5hr_up = drought05_up, drought_4.0hr_up =drough4_up,drought_0.5hr_down = drought05_down, drought_4hr_down =  drought4_down,
                   salt_0.5hr_up = salt05_up, salt_4.0hr_up = salt4_up, salt_0.5hr_down = salt05_down, salt_4hr_down = salt4_down, 
                   sporeSig = spore_sig, caulonemaSig = caul_sig, 
                   chloronemaSig = chlor_sig, gametophoreSig = gam_sig, rhizoidSig = rhiz_sig, 
                   sporeUp = spore_up, caulonemaUp = caul_up, chloronemaUp = chlor_up, gamUp = gam_up, rhizUp = rhiz_up)

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
posneg = c('pos','neg')
n = 0
for(i in 1:4){
  for(j in 1:2){
    n = n+1
    clustnames = c(clustnames, paste('cluster', i, posneg[j], sep = '_'))
    clustdat[[n]] = read.csv(paste('./sig_cluster_subsets/', posneg[j], '_c', i, '_6k.txt', sep = ''), header = FALSE)$V1
  }
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
  sig_compare[cl_map[[i]],i+length(names(rawsyn))] = TRUE
}
for(i in 1:length(maplist)){
  sig_compare[maplist[[i]], i+length(names(rawsyn)) + length(names(cl_map))] = TRUE
}

#it's all together! let's write it, and open a new script to work on comparing cols
write.csv('./outputs/12_significant_by_set_raw.txt', row.names = FALSE, x = sig_compare, quote = FALSE)
