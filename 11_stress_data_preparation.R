#prepare data from rpkm stress data file

#set working directory
setwd('/Users/Josophila/Desktop/SEQ_data/10x_version3_genome_reanalysis/')

stressFC_start = read.csv('../moss_omics_datasets/stress-fold_change_data.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
stressFC_one = stressFC_start[,-c(10, 11)]
stressFC_one[is.na(stressFC_one)] = FALSE       #change NA values to FALSE because they muck up the works otherwise

stressList = c('ABA', 'Cold', 'Drought', 'Salt')     #will iterate over these names
tpts = c('.0.5.h', '.4.0.h')                         #will iterate over these timepoints


stress_ups = list()                                  #make an empty list here 
for(i in 1:length(stressList)){                      # This loop will take all gene names with log2fc above .1
  for(j in 1:length(tpts)){
    nam =paste(stressList[i],tpts[j], sep = '')
    namup = paste(stressList[i],tpts[j], 'UP', sep = '')
    stress_ups[[namup]] = stressFC_one[stressFC_one[,nam] >= 0.1,1]
  }
}

stress_downs = list()
for(i in 1:length(stressList)){
  for(j in 1:length(tpts)){
    nam = paste(stressList[i], tpts[j], sep = '')
    namup = paste(stressList[i], tpts[j], 'DOWN', sep = '')
    stress_downs[[namup]] = stressFC_one[stressFC_one[,nam] <= -0.1, 1]
  }
}

#compile ups and downs into a table
upanddown = c(stress_ups, stress_downs)

#write all of these to a file
for(i in 1:length(upanddown)){
  fnam = names(upanddown)[i]
  write.csv(upanddown[[i]], file = paste('./outputs/stress_data/',fnam,'.txt', sep = ''), row.names = FALSE)
}
