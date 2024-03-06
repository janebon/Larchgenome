setwd('G:/Мой диск/Bioinform/Lar Annot/Repeats')

#----------------------------- Unspecified part -------------------------------
data <- fread('Larch.essential.Unspecified.corrected.out')
names(data) <- c('start', 'end', 'name', 'group')

fu <- as.data.frame(table(data$name))


DNA <- c('Mariner', 'hAT', 'MuDR', 'EnSpm', 'CACTA', 'piggyBac', 'P', 'Merlin', 'Harbinger', 'Transib', 'Novosib', 
         'Helitron','Polinton', 'Kolobok', 'Crypton', 'Sola', 'Zator', 'Ginger1', 'Ginger2', 'Ginger', 
         'Academ', 'Zisupton', 'IS3EU', 'ISL2EU', 'Dada', 'DNA', 'Chapaev', 'Vandal', 'TIR')

LTR <- c('Gypsy', 'Copia', 'Gymny', 'BEL', 'DIRS', 'ERV1', 'ERV2', 'ERV3', 'ERV4', 'LTR', 'ERV', 
         'PtAngelina', 'PtAppalachian', 'PtBastrop', 'PtConagree', 'PtCumberland', 'PtOuachita', 
         'PtOzark', 'PtPineywoods', 'PtTalladega', 'Corky', 'IFG')

nonLTR <- c('SINE', 'CRE', 'NeSL', 'R4', 'R2', 'R1', 'L1', 'L2', 'RTE', "I", 
'Jockey', 'CR1', 'Rex1', 'RandI', 'Penelope', 'Tx1', 'RTEX', 'Crack', 
'Nimb', 'Proto1', 'Proto2', 'Proto', 'RTETP', 'Hero', 'Tad', 'Loa', 'Ingi', 
'Outcast', 'Daphne', 'Ambal', 'Vingi', 'Kiri', 'NonLTR', 'LINE', 'tRNA', 'Utopia')

InVir <- 'Caulimovirus'

tRNA <- 'tRNA'

data$group <- ifelse(data$name %in% DNA, 'DNA', data$group)
data$group <- ifelse(data$name %in% LTR, 'LTR', data$group)
data$group <- ifelse(data$name %in% nonLTR, 'nonLTR', data$group)
data$group <- ifelse(data$name %in% InVir, 'Integrated Virus', data$group)
data$group <- ifelse(data$name %in% tRNA, 'tRNA', data$group)
data$group <- ifelse(data$name == 'Unknown', 'Unclassified', data$group)
data$group <- ifelse(data$group == 'Unspecified', 'Other', data$group)

data$len <- data$end - data$start

write.table(data, 'Larch.repeats1', quote=FALSE, row.names=FALSE, sep="\t")



#-------------------------------- Spesified part   -----------------------------

data <- fread('Larch.essential.Specified.corrected.out', fill=TRUE, header=F, skip=3)
names(data) <- c('start', 'end', 'name')

fu <- as.data.frame(table(data$name))


DNA <- c('Mariner', 'hAT', 'MuDR', 'EnSpm', 'CACTA', 'piggyBac', 'P', 'Merlin', 'Harbinger', 'Transib', 'Novosib', 
         'Helitron','Polinton', 'Kolobok', 'Crypton', 'Sola', 'Zator', 'Ginger1', 'Ginger2', 'Ginger', 
         'Academ', 'Zisupton', 'IS3EU', 'ISL2EU', 'Dada', 'DNA', 'Chapaev', 'Vandal', 'TIR')

LTR <- c('Gypsy', 'Copia', 'Gymny', 'BEL', 'DIRS', 'ERV1', 'ERV2', 'ERV3', 'ERV4', 'LTR', 'ERV', 
         'PtAngelina', 'PtAppalachian', 'PtBastrop', 'PtConagree', 'PtCumberland', 'PtOuachita', 
         'PtOzark', 'PtPineywoods', 'PtTalladega', 'Corky', 'IFG')

nonLTR <- c('SINE', 'CRE', 'NeSL', 'R4', 'R2', 'R1', 'L1', 'L2', 'RTE', "I", 
'Jockey', 'CR1', 'Rex1', 'RandI', 'Penelope', 'Tx1', 'RTEX', 'Crack', 
'Nimb', 'Proto1', 'Proto2', 'Proto', 'RTETP', 'Hero', 'Tad', 'Loa', 'Ingi', 
'Outcast', 'Daphne', 'Ambal', 'Vingi', 'Kiri', 'NonLTR', 'LINE', 'tRNA', 'Utopia')

Other <- c('Other', 'TE', 'Centromeric', 'IID2-12_AT', 'Host')

data$group <- '-'
data$group <- ifelse(data$name %in% DNA, 'DNA', data$group)
data$group <- ifelse(data$name %in% LTR, 'LTR', data$group)
data$group <- ifelse(data$name %in% nonLTR, 'nonLTR', data$group)
data$group <- ifelse(data$name %in% Other, 'Other', data$group)
data$group <- ifelse(data$name == 'rRNA', 'rRNA', data$group)
data$group <- ifelse(data$name == 'Satellite', 'Satellite', data$group)
data$group <- ifelse(data$name == 'SSR', 'Satellite', data$group)
data$group <- ifelse(data$name == 'Simple_repeat', 'Simple repeat', data$group)
data$group <- ifelse(data$name == 'Retro', 'Retro', data$group)
data$group <- ifelse(data$name == 'Low_complexity', 'Low complexity', data$group)
data$group <- ifelse(data$name == 'Uncategorized', 'Unclassified', data$group)
data$group <- ifelse(data$name == 'Unknown', 'Unclassified', data$group)

data$len <- data$end - data$start

write.table(data, 'Larch.repeats2', quote=FALSE, row.names=FALSE, sep="\t")

	
#-------------------------------- Join   ---------------------------------------

rep1 <- fread('Larch.repeats1')
rep2 <- fread('Larch.repeats2')

rep <- rbind(rep1, rep2)

c <- as.data.frame(table(rep$group))
l <- aggregate(rep$len, by=list(rep$group), FUN=sum)
r <- cbind(c,l)
r[,3] <- NULL
names(r) <- c('Group', 'Number', 'Len')
r$percentL <- r$Len/12340000000*100
r$precentN <- r$Number / nrow(rep)

cc <- as.data.frame(table(rep$name))
ll <- aggregate(rep$len, by=list(rep$name), FUN=sum)
rr <- cbind(cc,ll)
rr[,3] <- NULL
names(rr) <- c('Group', 'Number', 'Len')
rr$percentL <- rr$Len/12340000000*100
rr$precentN <- rr$Number / nrow(rep)


write.table(r, 'Repeats.Groups.tbl', quote=FALSE, row.names=FALSE, sep="\t")
write.table(rr, 'Repeats.Names.tbl', quote=FALSE, row.names=FALSE, sep="\t")

#--------------------------------- Stats  --------------------------------------
options(scipen = 999)


data2 <- read.table('Repeats.Groups.tbl', header = TRUE, sep = '\t')   # groups
data2$percentL <- round(data2$percentL, 4)
data2$precentN <- round(data2$precentN, 4)


data <- read.table('Repeats.Names.tbl', header = TRUE)     # names
data <- data[order(-data$Number),]
data$percentL <- round(data$percentL, 2)
data$precentN <- round(data$precentN, 2)
data$fu <- ' '
data$fu <- ifelse(data$Group %in% DNA, 'DNA', data$fu)
data$fu <- ifelse(data$Group %in% LTR, 'LTR', data$fu)
data$fu <- ifelse(data$Group %in% nonLTR, 'nonLTR', data$fu)

names(data)[1] <- 'Clade'
names(data)[6] <- 'Group'


#  Make table for Plot

tbl <- subset(data, percentL > 0)
tbl <- tbl[c(-16,-17), ]
tbl[7,c(2,3,4,5)] <- tbl[7,c(-1,-6)] + tbl[20,c(-1,-6)]
tbl <- tbl[-20,]
tbl[3,c(2,3,4,5)] <- tbl[3,c(-1,-6)] + tbl[18,c(-1,-6)]
tbl <- tbl[-18,]
tbl[13,c(2,3,4,5)] <- tbl[13,c(-1,-6)] + tbl[2,c(-1,-6)]
tbl <- tbl[-2,]
tbl <- tbl[order(tbl$Group),]

write.table(tbl, 'Plot.data.txt', quote=FALSE, row.names=FALSE, sep="\t")



#  Make Plot  -------------------------------------- ---------------------------

library(ggthemes)
library(ggplot2)
library(ggsci)

tbl <- read.table('Plot.data.txt', sep = '\t', header=TRUE)


png(file='Repeats.png', width = 10, height = 5, units = 'in', res = 800)

# scale_fill_tableau(palette = "Tableau 20", name='Superfamily')
# cols <- colorRampPalette(pal_tron()(7))(17)
# cols <- colorRampPalette(pal_startrek()(7))(17)
# cols <- colorRampPalette(pal_locuszoom()(7))(17)
cols <- colorRampPalette(pal_lancet()(8))(17)

ggplot(tbl, aes(x=Group, y=Len/1000, fill=Clade))+
  geom_bar(stat="identity", color = 'black')+
  labs(x='Group/Subclass', y='Size (kbp)', tag='A')+
  theme_few()+
  theme(axis.text=element_text(size=13), 
        axis.title=element_text(size=15,face="bold"))+
  xlim( 'Uncategorized','Low complexity','Simple repeat','Other Retro','non-LTR','LTR','DNA')+
  coord_flip()+
  scale_fill_manual(values = cols)


dev.off()




#------------------------------- 


rep1 <- fread('Larch.repeats1')
rep2 <- fread('Larch.repeats2')
rep <- rbind(rep1, rep2)

rep <- subset(rep1, name=='Copia'|name=='Gypsy'|name=='LINE'|name=='L1'|name=='I')
rep$name <- ifelse(rep$name == 'LINE', 'LINE/L1', rep$name)
rep$name <- ifelse(rep$name == 'L1', 'LINE/L1', rep$name)

cols = c('#F1CE63', '#499894', '#E15759', '#FF9D9A')

png(file='Repeats_copynumber.png', width = 10, height = 5, units = 'in', res = 800)

ggplot(rep, aes(x=len, fill=name))+
  geom_bar(stat='count')+
  facet_wrap(~ name, scales = "free")+
  labs(x='TE Length, bp', y='Copy Number')+
  theme_bw()+
  scale_fill_manual(values = cols)+
  theme(legend.position="none")
  
dev.off()

