#Determine similarity vs dissimilarity of each cluster to each other dataset

#set working directory
setwd('/Users/Josophila/Desktop/SEQ_data/10x_version3_genome_reanalysis/')

#read individual cluster information
singleClusts = list(clust1_all = read.csv('./outputs/13_single_cluster1_all_comparisons.txt'),
                    clust2_all = read.csv('./outputs/13_single_cluster2_all_comparisons.txt'),
                    clust3_all = read.csv('./outputs/13_single_cluster3_all_comparisons.txt'),
                    clust4_all = read.csv('./outputs/13_single_cluster4_all_comparisons.txt'))


conditions = c('cold_0.5', 'cold_4.0', 'aba_0.5','aba_4.0','salt_0.5','salt_4.0',
              'drought_0.5','drought_4.0', 'gam','rhiz','chlor','caul','spore')

#loop is going to take all comparisons of a cluster to a given category and determine a similarity score. Similarity is going to be sum of each OBSv score in the interaction, where inverse correlations are negative. Independence score is just the sum of each of the cells, not taking into account the correlation. 
results = list()
for(j in 1:length(singleClusts)){
  df = singleClusts[[j]]
  independence = c()
  similarity = c()
for(i in 1:length(conditions)){
    subject = grep(pattern = conditions[i], x = df$Comparison)
    similarity = c(similarity, sum(df[subject, 'X.Obs...Exp..2.Exp'] * df[subject, 'Correlation'])/length(subject))
    independence = c(independence, 
                     sum(df[subject, 'X.Obs...Exp..2.Exp'])/length(subject))
  }
  results[[j]] = data.frame(conditions,nonrandom = independence, similarity)
}
#show each comparison with the calculated scores

for(i in 1:4){
  nam = paste('./outputs/14_similarity_scores_cluster',i,'.txt', sep = '')
  write.csv(file = nam, results[[i]], row.names = FALSE, quote = FALSE)
}
