#figure out how to analyze the dataset of gene and DE or not in each cluster/set
#as of 1/6/17 - got it all working and just need to change the clusters as they were extracted to distinguish UP in a cluster vs DOWN in a cluster

#set working directory
setwd('/Users/Josophila/Desktop/SEQ_data/10x_version3_genome_reanalysis/')

sigDat = read.csv('./outputs/12_significant_by_set_raw.txt', sep = ',', header = TRUE, stringsAsFactors = FALSE)

sigDatOnly = sigDat[,-c(1,2)]
#change NA to FALSEs
sigDatOnly[is.na(sigDatOnly)] = FALSE

#number significant in each column independently
numSig = colSums(sigDatOnly)

#probability that it is significant in two categories is P(A) * P(B); 
#then multiply by number of genes to get expected

#single probabilities
sprob = numSig/nrow(sigDat)

#pairwise probabilities
pwProb = combn(sprob, 2, prod)

#pairwise probabilities times the number of genes
exp_overlap = pwProb*nrow(sigDatOnly)


#pairwise sums of every column with every other
pwSums = combn(ncol(sigDatOnly), 2L, function(pwSums) rowSums(sigDatOnly[pwSums]))
colnames(pwSums) = combn(paste(names(sigDatOnly), ' ', sep = ''), 2L, paste, collapse = "")
pwSums = data.frame(pwSums)
#write this, since it contains all of the overlaps
#write.csv('./overlap_analysis_01-SUMS.txt', row.names = FALSE, x = pwSums)

#Count the number of overlapping hits, ie the '2' values
twoCount = c()
for(i in 1:ncol(pwSums)){
  twoCount = c(twoCount, sum(pwSums[,i] == 2))
}
names(twoCount) = colnames(pwSums)

#data frame with observed and expected values for each comparison
overlap_dat = data.frame(colnames(pwSums), expected = exp_overlap, detected = twoCount)
#take out 'sig' values from gene atlas; up and down are all we care about
removemap = grep(overlap_dat$colnames.pwSums., pattern = 'Sig')
overlap_dat = overlap_dat[-removemap,]

#keep only comparisons of clusters to everything
overlap_clusters = overlap_dat[grep(ignore.case = TRUE, x=overlap_dat[,1], pattern = 'clust'),]
row.names(overlap_clusters) = NULL

write.csv(file = './outputs/13_clusters_overlap.txt', overlap_clusters, row.names = FALSE, quote = FALSE)

#calculate something like chi squared for cells - ((Obs-Exp)^2)/Exp
ovClustDif = ((overlap_clusters$detected - overlap_clusters$expected)^2)/overlap_clusters$expected

#make frame including that difference measure
overlap_mes = data.frame(overlap_clusters, ovClustDif)

#score whether observed higher or lower than expected
overexp = which(overlap_mes$detected > overlap_mes$expected)
underexp = which(overlap_mes$detected < overlap_mes$expected)
overlap_mes[overexp,'high_vs_low'] = 'Higher Than Expected'
overlap_mes[underexp, 'high_vs_low'] = 'Lower Than Expected'

#Score whether positively correlated or inversely correlated (whether upregulated goes with upregulated vs up going with downregulated)
poslist = grep(pattern = 'pos', overlap_mes$colnames.pwSums.)
neglist = grep(pattern = 'neg', overlap_mes$colnames.pwSums.)
downlist = grep(pattern = c('down|Down'), overlap_mes$colnames.pwSums.)
uplist = grep(pattern = c('up|Up'), overlap_mes$colnames.pwSums.)

#find the intersection!
posdown = intersect(poslist, downlist)
posup = intersect(poslist, uplist)
negdown = intersect(neglist, downlist)
negup = intersect(neglist, uplist)

#Gives a plus if it is up and up regulated in each group, gives a - if it is down and downregulated in each group
overlap_mes[c(posdown, negup),'correlation'] = -1
overlap_mes[c(posup,negdown),'correlation'] = 1

#Get rid of cluster vs cluster comparisons, as they might mess it up
cVsC= grep(pattern = 'cluster_[1-4]_[a-z]{3}.cluster_[1-4]_[a-z]{3}', overlap_mes$colnames.pwSums., useBytes = TRUE)
overlaps_noCvC = overlap_mes[-cVsC,]

#change names
namnam = c('Comparison', 'Expected_Overlap', 'Observed_Overlap','(Obs - Exp)^2/Exp', 'Abover or Below', 'Correlation')
colnames(overlaps_noCvC) = namnam
colnames(overlap_mes) = namnam

#calculate diff score for overlap_dat (one with all comparisons)
allVsAllDif = ((overlap_dat$detected-overlap_dat$expected)^2)/overlap_dat$expected
allVsAllDat = data.frame(overlap_dat, allVsAllDif)
colnames(allVsAllDat) = namnam[1:4]

#get cluster by cluster overlap data
spot = list()
for (i in 1:4){
  nam = paste('^cluster_', i, sep = '')
  spot[[i]] = grep(pattern = nam, overlaps_noCvC$Comparison, useBytes = TRUE)
}

clustByClust = list(c1only = overlaps_noCvC[spot[[1]],], 
                    c2only = overlaps_noCvC[spot[[2]],],
                    c3only = overlaps_noCvC[spot[[3]],],
                    c4only = overlaps_noCvC[spot[[4]],])

#export each single cluster data, with only top 20 and with all in order
for(i in 1:4){
  bwah = clustByClust[[i]]
  bwah = bwah[order(bwah$'(Obs - Exp)^2/Exp'),]
  tops = head(bwah, 20)
  foob = paste('./outputs/13_single_cluster_',i,'top_20.txt', sep = '')
  ooby = paste('./outputs/13_single_cluster',i,'_all_comparisons.txt', sep = '')
  #obar = paste('c',i,'only', sep = '')
  #foobydoo = clustByClust[obar]
  write.csv(file = foob, x = tops, row.names = FALSE, quote = FALSE)
  write.csv(file = ooby, x = bwah, row.names = FALSE, quote =FALSE)
}

#write files
#1) All comparisons
write.csv(file = './outputs/13_overlap_all_comparisons.txt', allVsAllDat, row.names = FALSE, quote = FALSE)

#2) Clusters compared to everything
write.csv(file = './outputs/13_clusters_to_all.txt', overlap_mes, row.names = FALSE, quote = FALSE)

#3) cluster to cluster comparisons removed
write.csv(file = './outputs/13_overlap_no_clustVsClust.txt', overlaps_noCvC, row.names = FALSE, quote = FALSE)


