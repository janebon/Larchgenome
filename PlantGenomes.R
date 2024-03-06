library(rgbif)        # To lookup names in the GBIF backbone taxonomy
library(inborutils)   # To wrap GBIF API data
library(ggplot2)
library(ggthemes)
library(ggsci)
library(ggrepel)
library(scales)
library(patchwork)

#######################################################
setwd('G:/Мой диск/Bioinform/Dissertation/Доп.мат')

dt <- read.table('Plant_genomes_19.09.2323.txt', sep = '\t', col.names = 'Sp')
dt2 <- as.data.frame(dt[!duplicated(dt[,'Sp']),])
names(dt2)[1] <- 'Sp'


info <- gbif_species_name_match(df = dt2, name = 'Sp', gbif_terms = c("family", "order", 'class', 'phylum', 'kingdom')) 
as.data.frame(table(info$class))


gymnosperm <- subset(info, class=='Pinopsida' | class=='Cycadopsida' | class=='Ginkgoopsida' | class=='Gnetopsida' )
write.table(gymnosperm, 'Plant_genomes_19.09.2323.gymnosperm.txt', quote=FALSE, row.names=FALSE, sep="\t")





##################################################

setwd('G:/Мой диск/Bioinform/Dissertation/Доп.мат')

dt <- read.table('Genomes_stat.txt', sep = '\t', header = TRUE)

info <- gbif_species_name_match(df = dt, name = 'Species', gbif_terms = c("family", "order", 'class')) 
info <- info[,c("Clade", "class", "order", "family", "Species", "genome.Mb", "gene", "repeats")]

write.table(info, 'Genomes_stat.txt', quote=FALSE, row.names=FALSE, sep="\t")
#######
setwd('G:/Мой диск/Bioinform/Dissertation/Доп.мат')

info <- read.table('Genomes_stat.txt', sep = '\t', header = TRUE)


pal <- c(Cycadopsida = 'firebrick1', 
         Ginkgoopsida = 'darkorange',
         Gnetopsida = 'darkgoldenrod1',
         Liliopsida = 'royalblue',
         Magnoliopsida = 'mediumorchid',
         Pinopsida = 'springgreen4')


# Gene number
p1 <- ggplot(info, aes(x=genome.Mb, y=gene, color=class))+
  geom_point(alpha=0.6, size=6)+
  geom_label_repel(aes(label = Species, fontface="italic"),
                   label.size=NA,
                   label.padding = 0.1,
                   box.padding   = 0.5, 
                   point.padding = 0.5,
                   max.overlaps = 8, 
                   segment.color = 'grey50',
                   show.legend  = FALSE) +
  labs(x='Размер геномной сборки, млн п.н.', y='Число генов', color='')+
  scale_x_continuous(labels = label_number(big.mark = " "))+
  scale_y_continuous(labels = label_number(big.mark = " "))+
  scale_color_manual(values = pal)

# Repeat %
p2 <- ggplot(info, aes(x=genome.Mb, y=repeats, color=class))+
  geom_point(alpha=0.6, size=6)+
  geom_label_repel(aes(label = Species, fontface="italic"),
                   label.size=NA,
                   label.padding = 0.1,
                   box.padding   = 0.5, 
                   point.padding = 0.5,
                   max.overlaps = 7, 
                   segment.color = 'grey50',
                   show.legend  = FALSE) +
  labs(x='Размер геномной сборки, млн п.н.', y='Доля повторов, %', color='')+
  scale_x_continuous(labels = label_number(big.mark = " "))+
  scale_y_continuous(labels = label_number(big.mark = " "))+
  scale_color_manual(values = pal)



p1 / p2 + 
  plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = 'A') & 
  theme_gdocs() + theme(axis.text=element_text(size=12),
                        axis.title=element_text(size=15, face="bold", color="gray40"),
                        legend.text = element_text(size=13, color="gray60"),
                        legend.position = "bottom", legend.direction = "horizontal",
                        axis.line.x.bottom = element_line(color="gray60", size = 1),
                        axis.line.y.left = element_line(color="gray60", size = 1),
                        plot.background = element_rect(color = "white"))
