#set working directory to SEQ/10x reanalysis.  In directory with DE csvs.
setwd('/Users/Josophila/Desktop/SEQ_data/10x_version3_genome_reanalysis/')

#import datasets
de3k = read.csv('./G04-3k_differential_expression.csv', header = TRUE, sep = ',', stringsAsFactors = FALSE)
de6k = read.csv('./G04-6K_differential_expression.csv', header=TRUE, sep=',', stringsAsFactors = FALSE)

head(de3k)
summary(de3k)

#generate matrices that will have all of the info we want summarized
Scoop3k = matrix(ncol = 5, nrow = 5); 
Scoop3k[1,1] = 'info gransden 3k cells'
Scoop6k = matrix(ncol = 5, nrow = 5); 
Scoop6k[1,1] = 'info gransden 6k cells'
cnames = c('cluster 1', 'cluster 2', 'cluster 3', 'cluster 4')
rnames = c('number of genes', 'number of DE genes', 'number up', 'number down')
Scoop3k[1,2:5] = cnames ; Scoop3k[2:5,1] = rnames
Scoop6k[1,2:5] = cnames ; Scoop6k[2:5,1] = rnames

#find the number of genes detected in each cluster
c1_3detected = sum(de3k$Cluster.1.Mean.UMI.Counts > 0)
c2_3detected = sum(de3k$Cluster.2.Mean.UMI.Counts > 0)
c3_3detected = sum(de3k$Cluster.3.Mean.UMI.Counts > 0)
c4_3detected = sum(de3k$Cluster.4.Mean.UMI.Counts > 0)

c1_6detected = sum(de6k$Cluster.1.Mean.UMI.Counts > 0)
c2_6detected = sum(de6k$Cluster.2.Mean.UMI.Counts > 0)
c3_6detected = sum(de6k$Cluster.3.Mean.UMI.Counts > 0)
c4_6detected = sum(de6k$Cluster.4.Mean.UMI.Counts > 0)

#load into vector, then matrix
detected3 = c(c1_3detected, c2_3detected, c3_3detected, c4_3detected)
detected6 = c(c1_6detected, c2_6detected, c3_6detected, c4_6detected)

Scoop3k[2,2:5] = detected3
Scoop6k[2,2:5] = detected6

#find number of differentially expressed gene (p< 0.01) in each cluster
c1_3DE = sum(de3k$Cluster.1.Adjusted.p.value < 0.01)
c2_3DE = sum(de3k$Cluster.2.Adjusted.p.value < 0.01)
c3_3DE = sum(de3k$Cluster.3.Adjusted.p.value < 0.01)
c4_3DE = sum(de3k$Cluster.4.Adjusted.p.value < 0.01)

c1_6DE = sum(de6k$Cluster.1.Adjusted.p.value < 0.01)
c2_6DE = sum(de6k$Cluster.2.Adjusted.p.value < 0.01)
c3_6DE = sum(de6k$Cluster.3.Adjusted.p.value < 0.01)
c4_6DE = sum(de6k$Cluster.4.Adjusted.p.value < 0.01)

de3_clusts = c(c1_3DE, c2_3DE, c3_3DE, c4_3DE)
de6_clusts = c(c1_6DE, c2_6DE, c3_6DE, c4_6DE)

Scoop3k[3,2:5] = de3_clusts
Scoop6k[3,2:5] = de6_clusts

#Make subsets of significant values
siglist3c1 = subset(de3k, de3k$Cluster.1.Adjusted.p.value < 0.01)
siglist3c2 = subset(de3k, de3k$Cluster.2.Adjusted.p.value < 0.01)
siglist3c3 = subset(de3k, de3k$Cluster.3.Adjusted.p.value < 0.01)
siglist3c4 = subset(de3k, de3k$Cluster.4.Adjusted.p.value < 0.01)

siglist6c1 = subset(de6k, de6k$Cluster.1.Adjusted.p.value < 0.01)
siglist6c2 = subset(de6k, de6k$Cluster.2.Adjusted.p.value < 0.01)
siglist6c3 = subset(de6k, de6k$Cluster.3.Adjusted.p.value < 0.01)
siglist6c4 = subset(de6k, de6k$Cluster.4.Adjusted.p.value < 0.01)


#Detect number up
positive3 = c(sum(siglist3c1$Cluster.1.Log2.fold.change > 0), 
              sum(siglist3c2$Cluster.2.Log2.fold.change > 0),
              sum(siglist3c3$Cluster.3.Log2.fold.change > 0),
              sum(siglist3c4$Cluster.4.Log2.fold.change > 0))
positive6 = c(sum(siglist6c1$Cluster.1.Log2.fold.change > 0),
              sum(siglist6c2$Cluster.2.Log2.fold.change > 0),
              sum(siglist6c3$Cluster.3.Log2.fold.change > 0),
              sum(siglist6c4$Cluster.4.Log2.fold.change > 0))

negative3 = c(sum(siglist3c1$Cluster.1.Log2.fold.change < 0), 
              sum(siglist3c2$Cluster.2.Log2.fold.change < 0),
              sum(siglist3c3$Cluster.3.Log2.fold.change < 0),
              sum(siglist3c4$Cluster.4.Log2.fold.change < 0))
negative6 = c(sum(siglist6c1$Cluster.1.Log2.fold.change < 0),
              sum(siglist6c2$Cluster.2.Log2.fold.change < 0),
              sum(siglist6c3$Cluster.3.Log2.fold.change < 0),
              sum(siglist6c4$Cluster.4.Log2.fold.change < 0))

Scoop3k[4,2:5] = positive3
Scoop3k[5,2:5] = negative3
Scoop6k[4,2:5] = positive6
Scoop6k[5,2:5] = negative6


write.csv(Scoop3k, file = './outputs/00_3ktable.csv', row.names = FALSE)
write.csv(Scoop6k, file = './outputs/00_6ktable.csv', row.names = FALSE)
