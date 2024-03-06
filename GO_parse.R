library(stringr)
library(tidyr)

setwd('G:/Мой диск/Bioinform/Lar Annot/B2G/parsing_GO/raw_files')

LIST1=list.files(pattern = '*.txt')
LIST1 = gsub('_gene_sets', '', LIST1)
LIST1 = unique(LIST1)

LIST2=list.files(pattern = '*gene_sets.txt')



for(i in 1:length(LIST1)){
   file=LIST1[i]
   sets=LIST2[i]

  # Get gene counts for GO IDs
  geneset <- read.table(sets, fill = TRUE, sep = '\t', header = TRUE)
  names(geneset)[1] <- 'GO.IDs'
  
  geneset$count <- str_count(geneset$Sequences, "\\S+")                         # count number of genes listed for this GO
  geneset$Sequences <- NULL

  # Get descriptions for GO IDs
  gos <- read.table(file, sep = '\t', fill = TRUE, header = TRUE, comment.char = '', quote='')
  gos <- subset(gos, Tags == '[BLASTED, MAPPED, ANNOTATED]')            # select only annotated genes
  gos <- subset(gos, GO.IDs != '')                                      # omit genes with emptly GO IDs
  total <- length(unique(gos$SeqName))                                  # Get total number of genes annotated with GO (for percent)
  gos <- gos[, c('GO.IDs', 'GO.Names')]                                 # select columns with ID and Description info
  
  gos2 <- separate_rows(gos, GO.IDs, sep = '; ')                        # split IDs 
  gos3 <- separate_rows(gos, GO.Names, sep = '; ')                      # split Descriptions
  
  descr <- cbind(gos2[,1], gos3[,2])                                    # join IDs and Descriptions
  descr <- subset(descr, GO.IDs != '')                                  # omit empty rows
  descr <- unique(descr)                                                # remove duplicates (there are many)
  rm(gos2,gos3, gos)
  
  descr$Ontology <- sub(":.*", "", descr$GO.Names)                      # add Ontology column
  descr$GO.Names <- gsub("P:|F:|C:", "", descr$GO.Names, perl=TRUE)     # remove Ontology indicator
  descr$GO.IDs <- gsub("P:|F:|C:", "", descr$GO.IDs, perl=TRUE)

  # Merge Counts and Descriptions
  
  dt <- merge(geneset, descr)
  dt <- dt[, c('GO.IDs', 'GO.Names', 'Ontology', 'count')]
  dt$perc <- round(dt$count / total * 100, 3)
  

  write.table(dt, paste0('G:/Мой диск/Bioinform/Lar Annot/B2G/parsing_GO/processed_files/', file, '.processed'), quote=FALSE, row.names=FALSE, sep="\t")
  write.table(paste(file, total), 'G:/Мой диск/Bioinform/Lar Annot/B2G/parsing_GO/processed_files/total', quote=FALSE, row.names=FALSE, sep="\t", append=T)

  }





