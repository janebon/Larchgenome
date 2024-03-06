setwd('D:/')

library(ggplot2)
library(ggsci)
library(RColorBrewer)




#### Select categories --------

L <- read.table('G:/Мой диск/Bioinform/Lar Annot/B2G/parsing_GO/processed_files/Lar_noIP.txt.processed', header = TRUE, sep = "\t", quote='')
geneset <-  read.table('G:/Мой диск/Bioinform/Lar Annot/B2G/parsing_GO/raw_files/Lar_noIP_gene_sets.txt', header = TRUE, sep = "\t", quote='')
names(geneset)[1] <- 'GO.IDs'

CC <- subset(L, Ontology == 'C')

cell <- CC[grep("cell", CC$GO.Names), ] 
fu <- merge(cell, geneset)
write.table(fu$Sequences , 'D:/genes_cell.IDs.txt', quote=FALSE, row.names=FALSE, sep="\t")

mito <- CC[grep("mitochondri", CC$GO.Names), ] 
fu <- merge(mito, geneset)
write.table(fu$Sequences , 'D:/genes_mito.IDs.txt', quote=FALSE, row.names=FALSE, sep="\t")

chloro <- CC[grep("chloroplast|plastid", CC$GO.Names), ] 
fu <- merge(chloro, geneset)
write.table(fu$Sequences , 'D:/genes_chloro.IDs.txt', quote=FALSE, row.names=FALSE, sep="\t")

nucl <- CC[grep("nucleus|nuclear", CC$GO.Names), ] 
fu <- merge(nucl, geneset)
write.table(fu$Sequences , 'D:/genes_nucl.IDs.txt', quote=FALSE, row.names=FALSE, sep="\t")

membr <- CC[grep("membrane", CC$GO.Names), ] 
fu <- merge(membr, geneset)
write.table(fu$Sequences , 'D:/genes_membr.IDs.txt', quote=FALSE, row.names=FALSE, sep="\t")

other <- sum(CC$count) - (sum(cell$count) + sum(mito$count) + sum(chloro$count) + sum(nucl$count) + sum(membr$count))

tbl_CC <- data.frame(ontology = 'CC',
                  category = c('cell', 'mito', 'chloro', 'nucl', 'membrane', 'other'),
                  count = c(sum(cell$count), sum(mito$count), sum(chloro$count), sum(nucl$count), sum(membr$count), other))




MF <- subset(L, Ontology == 'F')

bind <- MF[grep("binding", MF$GO.Names), ] 
fu <- merge(bind, geneset)
write.table(fu$Sequences , 'D:/genes_bind.IDs.txt', quote=FALSE, row.names=FALSE, sep="\t")

transcrip <- MF[grep("transcription", MF$GO.Names), ] 
fu <- merge(transcrip, geneset)
write.table(fu$Sequences , 'D:/genes_transcription.IDs.txt', quote=FALSE, row.names=FALSE, sep="\t")

oxi <- MF[grep("reductase", MF$GO.Names), ] 
fu <- merge(oxi, geneset)
write.table(fu$Sequences , 'D:/genes_oxidoreduct.IDs.txt', quote=FALSE, row.names=FALSE, sep="\t")

trans <- MF[grep("transferase|kinase|transaminase", MF$GO.Names), ] 
fu <- merge(trans, geneset)
write.table(fu$Sequences , 'D:/genes_transferase.IDs.txt', quote=FALSE, row.names=FALSE, sep="\t")

hydr <- MF[grep("nuclease|hydrolase|esterase|phosphatase|phospholipase|ribonuclease|Chitinase|helicase|atpase", MF$GO.Names), ] 
fu <- merge(hydr, geneset)
write.table(fu$Sequences , 'D:/genes_hydrolase.IDs.txt', quote=FALSE, row.names=FALSE, sep="\t")

transport <- MF[grep("transporter", MF$GO.Names), ] 
fu <- merge(transport, geneset)
write.table(fu$Sequences , 'D:/genes_transport_MF.IDs.txt', quote=FALSE, row.names=FALSE, sep="\t")

other = sum(MF$count) - (sum(bind$count)+ sum(transcrip$count)+ sum(oxi$count)+ sum(trans$count)+ sum(hydr$count)+ sum(transport$count))


tbl_MF <- data.frame(ontology = 'MF',
                  category = c('binding', 'transcription', 'oxidoreductase', 'transferase', 'hydrolase', 'transporter', 'other'),
                  count = c(sum(bind$count), sum(transcrip$count), sum(oxi$count), 
                            sum(trans$count), sum(hydr$count), sum(transport$count), other))



BP <- subset(L, Ontology == 'P')

cellp <- BP[grep("cellular", BP$GO.Names), ] 
fu <- merge(cellp, geneset)
write.table(fu$Sequences , 'D:/genes_cellularproc.IDs.txt', quote=FALSE, row.names=FALSE, sep="\t")

dev <- BP[grep("development", BP$GO.Names), ] 
fu <- merge(dev, geneset)
write.table(fu$Sequences , 'D:/genes_dev.IDs.txt', quote=FALSE, row.names=FALSE, sep="\t")

metabol <- BP[grep("metabol", BP$GO.Names), ] 
fu <- merge(metabol, geneset)
write.table(fu$Sequences , 'D:/genes_metabl.IDs.txt', quote=FALSE, row.names=FALSE, sep="\t")

regul <- BP[grep("regulation", BP$GO.Names), ] 
fu <- merge(regul, geneset)
write.table(fu$Sequences , 'D:/genes_regul.IDs.txt', quote=FALSE, row.names=FALSE, sep="\t")

response <- BP[grep("response", BP$GO.Names), ] 
fu <- merge(response, geneset)
write.table(fu$Sequences , 'D:/genes_response.IDs.txt', quote=FALSE, row.names=FALSE, sep="\t")

signal <- BP[grep("signal", BP$GO.Names), ] 
fu <- merge(signal, geneset)
write.table(fu$Sequences , 'D:/genes_signal.IDs.txt', quote=FALSE, row.names=FALSE, sep="\t")

biosynt <- BP[grep("biosynthe", BP$GO.Names), ] 
fu <- merge(biosynt, geneset)
write.table(fu$Sequences , 'D:/genes_biosynth.IDs.txt', quote=FALSE, row.names=FALSE, sep="\t")

catab <- BP[grep("catabolic", BP$GO.Names), ] 
fu <- merge(catab, geneset)
write.table(fu$Sequences , 'D:/genes_catabol.IDs.txt', quote=FALSE, row.names=FALSE, sep="\t")

transport <- BP[grep("transport", BP$GO.Names), ] 
fu <- merge(transport, geneset)
write.table(fu$Sequences , 'D:/genes_transport_BP.IDs.txt', quote=FALSE, row.names=FALSE, sep="\t")

other = sum(BP$count) - (sum(cellp$count)+ sum(dev$count)+ sum(metabol$count)+ sum(biosynt$count)+ sum(catab$count) +
                           sum(regul$count)+ sum(response$count)+ sum(signal$count)+ sum(transport$count))


tbl_BP <- data.frame(ontology = 'BP',
                     category = c('cellp', 'dev', 'metabol', 'regul', 'response', 'signal', 'biosynthetic', 'catabolic', 'transport', 'other'),
                     count = c(sum(cellp$count), sum(dev$count), sum(metabol$count), sum(regul$count), 
                               sum(response$count), sum(signal$count), sum(biosynt$count), sum(catab$count), sum(transport$count),other))


#### Make categories plot for Larch -------

tbl <- rbind(tbl_CC, tbl_MF, tbl_BP)


tbl_BP$category <- as.factor(tbl_BP$category)
tbl_BP$category <- factor(tbl_BP$category, levels = c('cellp', 'dev', 'metabol', 'regul', 'response', 'signal', 'biosynthetic', 'catabolic', 'transport', 'other'))


tbl_MF$category <- as.factor(tbl_MF$category)
tbl_MF$category <- factor(tbl_MF$category, levels = c('binding', 'transcription', 'oxidoreductase', 'transferase', 'hydrolase', 'transporter', 'other'))


tbl_CC$category <- as.factor(tbl_CC$category)
tbl_CC$category <- factor(tbl_CC$category, levels = c('cell', 'mito', 'chloro', 'nucl', 'membrane', 'other'))


cols <- colorRampPalette(brewer.pal(9, 'Blues'))(10) 

pBP <- ggplot(tbl_BP, aes(x=ontology, y=count, fill=category))+
  geom_bar(position="stack", stat="identity", color='black')+
  coord_flip()+
  geom_text(data=tbl_BP, aes(label = category), position = position_stack(vjust = 0.5), check_overlap=TRUE)+
  theme_classic()+
  scale_fill_manual(values = cols)


pMF <- ggplot(tbl_MF, aes(x=ontology, y=count, fill=category))+
  geom_bar(position="stack", stat="identity", color='black')+
  coord_flip()+
  geom_text(data=tbl_MF, aes(label = category), position = position_stack(vjust = 0.5), check_overlap=TRUE)+
  theme_classic()+
  scale_fill_brewer(palette='Purples')


pCC <- ggplot(tbl_CC, aes(x=ontology, y=count, fill=category))+
  geom_bar(position="stack", stat="identity", color='black')+
  coord_flip()+
  geom_text(data=tbl_CC, aes(label = category), position = position_stack(vjust = 0.5), check_overlap=TRUE)+
  theme_classic()+
  scale_fill_brewer(palette='Oranges')



library(patchwork)


svg("Rplots.svg" ,  width = 12, height = 7)

pBP / pMF / pCC

dev.off()



write.table(tbl, 'G:/Мой диск/Bioinform/Lar Annot/B2G/parsing_GO/Fun_summary.txt', quote=FALSE, row.names=FALSE, sep="\t")






######## parse BLAST results --------------
setwd('D:/larch')


## Count number of unique genes per category  ------

setwd('D:/larch/Arabid_compare/category_gene_lists')
FILES=list.files(pattern = 'txt.2.3')

total <- data.frame(cat = character(),
                  total = character())

for(i in 1:length(FILES)){
  f <- read.table(FILES[i], header = FALSE, sep = "\t", quote='')
  
  g <- unique(f)
  
  fu <- data.frame(cat = paste(gsub("\\.IDs.*$","", FILES[i])),
                      total = paste(nrow(g)))
  
  total <- merge(total, fu, all=TRUE)
  
  #write.table(g, paste(FILES[i], '.4', sep = ''), quote=FALSE, row.names=FALSE, sep="\t")  # save unique genes to get their fastas
  
}



## Count number of mapped genes per category -- pident > 50 & qcovhsp > 50  ----

setwd('D:/larchD/Arabid_compare')
FILES=list.files(pattern = '*.psl')

mapped <- data.frame(cat = character(),
                    mapped = character())

for(i in 1:length(FILES)){

  f <- read.table(FILES[i], header = FALSE, sep = "\t", quote='', comment.char = "")
  names(f) <- c('qseqid', 'stitle', 'qlen', 'slen', 'length', 'qcovhsp', 'evalue', 'pident')
  
  filtr <- subset(f, pident > 50 & qcovhsp > 50)
  un <- filtr[!duplicated(filtr[,'qseqid']),]
  
  fu <- data.frame(cat = paste(gsub("\\.IDs.*$","", FILES[i])),
                   mapped = paste(nrow(un)))
  
  mapped <- merge(mapped, fu, all = TRUE)
}

rm(fu, g, f, filtr, un)

# make table for plotting
dt <- merge(total, mapped)
dt$total <- as.numeric(dt$total)
dt$mapped <- as.numeric(dt$mapped)
dt$perc <- round(dt$mapped / dt$total *100)
dt$ontology <- as.factor(dt$ontology)

write.table(dt, 'Arabid_dt.txt', quote=FALSE, row.names=FALSE, sep="\t")




## Make plot ----

dt <- read.table('D:/larch/Arabid_compare/Arabid_dt.txt', header = TRUE, sep = '\t')

fix(dt)   # add ontology column

dt_BP <- subset(dt, ontology == 'BP')
dt_MF <- subset(dt, ontology == 'MF')
dt_CC <- subset(dt, ontology == 'CC')

dt_BP$cat <- as.factor(dt_BP$cat)
dt_MF$cat <- as.factor(dt_MF$cat)
dt_CC$cat <- as.factor(dt_CC$cat)

dt_BP$cat <- factor(dt_BP$cat, levels = c('Cellular process', 'Development', 'Metabolism', 'Regulation', 'Response to stimuli', 'Signaling', 'Biosynthesis', 'Catabolism', 'Transport'))
dt_MF$cat <- factor(dt_MF$cat, levels = c('Binding activity', 'Transcription activity', 'Oxidoreductase activity', 'Transferase activity', 'Hydrolase activity', 'Transporter activity'))
dt_CC$cat <- factor(dt_CC$cat, levels = c('Cell', 'Mitochondrion', 'Chloroplast', 'Nucleus', 'Membrane'))


library(RColorBrewer)
library(ggplot2)
library(patchwork)

colsBP =  colorRampPalette(c('ghostwhite','#466f9d'))(10)
colsMF =  colorRampPalette(c('mintcream','seagreen'))(6)
colsCC =  colorRampPalette(c('mistyrose','#ed444a'))(5)

p1 <- ggplot(dt_BP, aes(x = cat, y = perc, fill = cat))+
  geom_bar(stat = 'identity', alpha = 1, color = 'black')+
  scale_fill_manual(values = colsBP)+
  labs(y = 'Percent of mapped genes')+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        legend.position="none")

p2 <- ggplot(dt_MF, aes(x = cat, y = perc, fill = cat))+
  geom_bar(stat = 'identity', alpha = 1, color = 'black')+
  scale_fill_manual(values = colsMF)+
  labs(y = 'Percent of mapped genes')+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position="none")

p3 <- ggplot(dt_CC, aes(x = cat, y = perc, fill = cat))+
  geom_bar(stat = 'identity', alpha = 1, color = 'black')+
  scale_fill_manual(values = colsCC)+
  labs(y = 'Percent of mapped genes')+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position="none")



svg("Categories_Arabid.svg" ,  width = 12, height = 4)

p1 + p2 + p3 & theme(axis.text.x = element_text(size = 11, angle = 40, hjust = 1))

dev.off()


png("Categories_Arabid.png" ,  width = 12, height = 4, units = 'in', res = 600)

p1 + p2 + p3 & theme(axis.text.x = element_text(size = 11, angle = 40, hjust = 1))

dev.off()