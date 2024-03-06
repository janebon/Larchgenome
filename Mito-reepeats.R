setwd('C:/Users/Пользователь/Google Диск/Bioinform/Dissertation/Доп.мат/mito')

setwd('G:/Мой диск/Bioinform/Dissertation/mito')


library(data.table)
library(ggplot2)
library(ggsci)
library(ggthemes)

library(tidyverse)

rps <- fread('7_scaffolds.fa-2.out', fill=T)  # Load data
rps[, c(1,2,3,4,5,8,9,12,13,14,15,16)] <- NULL  # remove extra columns
rps$len <- rps$end - rps$begin  # add Lengtfh column
rps[,c(1,2)] <- NULL
colnames(rps) <- c("Repeat", "Class", 'Len')
rps$Class <- as.factor(rps$Class)

rps$Family <- rps$Class  # add Family column
rps$Class <- gsub("/.*$", "", as.matrix(rps$Class), perl=TRUE) # remove Family info from Class column
rps$Family <- gsub("non-LTR/RTEX", "non-LTR", as.matrix(rps$Family), perl=TRUE)
rps$Family <- gsub("Host", "Unknown", as.matrix(rps$Family), perl=TRUE)
rps$Class <- gsub("Host", "Unknown", as.matrix(rps$Class), perl=TRUE)
rps$Family <- gsub("DNA/.*$", "DNA", as.matrix(rps$Family), perl=TRUE) # remove DNA families
rps$Family <- gsub("LTR/DIRS", "LTR", as.matrix(rps$Family), perl=TRUE) # 
rps$Family <- gsub("non-LTR/I", "non-LTR", as.matrix(rps$Family), perl=TRUE) # 
rps$Family <- gsub("non-LTR/Penelope", "non-LTR", as.matrix(rps$Family), perl=TRUE) # 
rps$Class <- gsub("Retro_unclassified", "Other Retro", as.matrix(rps$Class), perl=TRUE)
rps$Family <- gsub("Retro_unclassified", "Other Retro", as.matrix(rps$Family), perl=TRUE)
rps$Family <- ifelse(rps$Repeat == "Gypsy", "LTR/Gypsy", rps$Family)
rps$Class <- ifelse(rps$Repeat == "Gypsy", "LTR", rps$Class)
rps$Family <- ifelse(rps$Repeat == "Copia", "LTR/Copia", rps$Family)
rps$Class <- ifelse(rps$Repeat == "Copia", "LTR", rps$Class)

rps2 <- rps %>% 
  group_by(Class, Family) %>% 
  summarise(Len = sum(Len))

rps2$percnt <- round(rps2$Len * 100 / 11690000, digits = 1)

write.table(rps2, 'Repeats_summary.txt', quote=FALSE, row.names=FALSE, sep="\t")



# Plot repeats Length by Class and color by Family
cols <- colorRampPalette(pal_startrek()(7))(12)
cols <- colorRampPalette(pal_npg()(8))(17)
cols <- colorRampPalette(pal_lancet()(8))(12)


rps2 <- read.table('Repeats_summary.txt', sep='\t', header = TRUE)

options(OutDec= ",")

svg("Mito_Repeats-2.svg", width = 10, height = 5)


ggplot(rps2, aes(x=Class, y=Len/1000, fill=Family))+
  geom_bar(stat="identity", color = 'black')+
  geom_text(aes(label = paste0(percnt, "%")), position = position_stack(vjust = 0.5), size = 4.5) +
  labs(x='Group/Subclass', y='Size (kbp)')+
  theme_few()+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15,face="bold"))+
  xlim( 'Unknown', 'rRNA', 'Low_complexity','Simple_repeat', 'Other Retro','non-LTR','LTR','DNA')+
  coord_flip()+
  scale_fill_tableau(palette = "Red-Blue-Brown", breaks=c('DNA' ,'LTR', "LTR/Copia", "LTR/Gypsy", 'non-LTR', "non-LTR/LINE", "non-LTR/SINE",
  'Other Retro','Simple_repeat', 'Low_complexity', 'rRNA','Unknown'))


dev.off()



