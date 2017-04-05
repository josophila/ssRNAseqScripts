#running gotime function on each file in directory
#seems like it could be best to run with full set of DE genes from each cluster

#set working directory
setwd('/Users/Josophila/Desktop/SEQ_data/10x_version3_genome_reanalysis/')

source('./scriptify-topgo.R')

#gotime(geneFile = './sig_cluster_totals/c1_3k_allsig.txt', outputFile = './TEST', descript = 'test')
#Do all three ontologies 
onto = c("BP", "MF", "CC")

for(k in 1:3){
  ont = onto[k]
  for(i in 1:4){
    for(j in 1:2){
      #iterate over loop to call different file names
     foo = paste('./sig_cluster_totals/c', i,'_',3*j,'k_allsig.txt', sep ='')
     gotime(geneFile = foo, 
           outputFile = paste('./top_GO_20_outputs/', 'topgo_20term_C', 
                             i,'_', ont, '_', 3*j, 'k.txt', sep = ''),
            descript = 'cluster1, 3k, top 100', set_type = ont, numterms = 20)
   }
  }
}
