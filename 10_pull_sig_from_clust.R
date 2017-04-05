#export_top sigs of each cluster
#Working in directory with DE 6k and DE 3k tables

#import tables
sixk = read.csv('./G04-6K_differential_expression.csv', sep = ',', header = TRUE, stringsAsFactors = FALSE)
threek = read.csv('./G04-3k_differential_expression.csv', sep = ',', header = TRUE, stringsAsFactors = FALSE)

#THIS IS REALLY HEAD - it ordered things in the opposite manner to what I expected. 
#I wound up taking the LEAST significant genes instead of the most. sigClustAll is better anyway for GO

sigclustTail = function(cluster, listnum = 50, frame){
  foo = paste('Cluster.', cluster ,'.Adjusted.p.value', sep = '')
  head(frame[order(frame[,foo]),], n = listnum)$Gene.ID  #don't worry about it... :(
}

sigClustAll = function(cluster, pval = 0.01, frame){
  foo = paste('Cluster.', cluster,'.Adjusted.p.value', sep = '')
  subset(frame, frame[,foo] < pval)$Gene.ID
}

sigClustpos = function(cluster, pval = 0.01, frame){
  foo = paste('Cluster.', cluster,'.Adjusted.p.value', sep = '')
  bar = paste('Cluster.', cluster,'.Log2.fold.change', sep = '')
  subset(frame, frame[,foo] < pval & frame[,bar] > 1)$Gene.ID
}
sigClustneg = function(cluster, pval = 0.01, frame){
  foo = paste('Cluster.', cluster,'.Adjusted.p.value', sep = '')
  bar = paste('Cluster.', cluster,'.Log2.fold.change', sep = '')
  subset(frame, frame[,foo] < pval & frame[,bar] < -1)$Gene.ID
}

# c1_3k_50 = sigclustTail(1, listnum = 50, frame = threek)
# write(c1_3k_50, './cluster1_3k_50sig.txt')
# #pull out top 50 - 200 genes from each cluster of 3k data
# for(i in 1:4){
#   for(j in 1:4){
#     blah = sigclustTail(cluster = i, listnum = 50*j, frame = threek)
#     write(blah, file = paste('./c',i,'k-3-',50*j, '.txt', sep = ''))
#   }
# }

#pull out top 50 - 200 genes from each cluster of 6k data
for(i in 1:4){
  for(j in 1:4){
    blah = sigclustTail(cluster = i, listnum = 50*j , frame = sixk)
    write(blah, file = paste('./c',i,'k-6-',50*j, '.txt',sep=''))
  }
}

#pull out all sig genes from each cluster of 3k data
for(i in 1:4){
  blah = sigClustAll(cluster = i, frame = threek)
  write(blah, file = paste('./c',i,'_3k_allsig.txt', sep =''))
}

#pull out all sig genes from each cluster of 6k data
for(i in 1:4){
  blah = sigClustAll(cluster = i, frame = sixk)
  write(blah, file = paste('./c',i,'_6k_allsig.txt',sep=''))
}

#pull out all sig and upregulated genes from each cluster of 6k data
for(i in 1:4){
  blah = sigClustpos(cluster = i, frame = sixk)
  write(blah, file = paste('./sig_cluster_subsets/pos_c', i,'_6k.txt', sep=''))
}
for(i in 1:4){
  blah = sigClustneg(cluster = i, frame = sixk)
  write(blah, file = paste('./sig_cluster_subsets/neg_c', i,'_6k.txt', sep=''))
}
