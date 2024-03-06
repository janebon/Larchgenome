
####### Make total gene counts for BP/MF/CC in each species -------

setwd('D:/GAF')

FILES=list.files(pattern = '.gaf') 
res <- data.frame(sp='', 
                  BP='',
                  CC='',
                  MF='',
                  stringsAsFactors=FALSE)
total <- data.frame(sp=character(), 
                    BP=character(),
                    CC=character(),
                    MF=character(),
                    stringsAsFactors=FALSE)

for (i in 1:length(FILES)){
  
  f = FILES[i]
  data <- read.table(f, skip = 4, sep = '\t')
  data <- data[, c(2, 5, 9, 10)] ; names(data) <- c('gene', 'term', 'domain', 'def')
  
  BP <- subset(data, domain == 'P')
  CC <- subset(data, domain == 'C')
  MF <- subset(data, domain == 'F')
  
  res$sp = gsub("_(.*?)$","", f)
  res$BP = length(unique(BP$gene))
  res$CC = length(unique(CC$gene))
  res$MF = length(unique(MF$gene))
  
  total <- merge(res, total, all=TRUE)
  
}  

total$total <- rowSums( total[,2:4] )


setwd('G:/Мой диск/Bioinform/Lar Annot/B2G/parsing_GO/')
write.table(total, 'Total_per_domain.txt', quote=FALSE, row.names=FALSE, sep="\t")
