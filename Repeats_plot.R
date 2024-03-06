
setwd('G:/Мой диск/Bioinform/Dissertation/Доп.мат/Plots/repeats')


library(patchwork)
library(ggthemes)
library(ggplot2)
library(ggsci)

tbl <- read.table('Plot.data.txt', sep = '\t', header=TRUE)
data <- read.table('Plot.data.micro.txt', sep = '\t', header=TRUE)
data2 <- read.table('total_TRF_GMATO.txt', sep = '\t', header=TRUE)




# scale_fill_tableau(palette = "Tableau 20", name='Superfamily')
# cols <- colorRampPalette(pal_tron()(7))(17)
# cols <- colorRampPalette(pal_startrek()(7))(17)
# cols <- colorRampPalette(pal_locuszoom()(7))(17)
cols <- colorRampPalette(pal_lancet()(8))(17)


percnt <- tbl$percentL


p1 <- ggplot(tbl, aes(x=Group, y=Len/1000, fill=Clade))+
  geom_bar(stat="identity", color = 'black')+
  geom_text(aes(label = paste0(percnt, "%")), position = position_stack(vjust = 0.5), size = 4.5) +
  labs(x='Группа/Класс', y='Размер, тыс п.н.', tag='A')+
  theme_classic()+
  xlim( 'Uncategorized','Low complexity','Simple repeat','Other Retro','non-LTR','LTR','DNA')+
  coord_flip()+
  scale_fill_manual(values = cols)


p2 <- ggplot(data, aes(x=Motif, y=Density, fill=Species))+
  geom_bar(stat='identity',position="dodge", color = 'black')+
  xlim('di', 'tri', 'tetra', 'penta', 'hexa', 'hepta', 'octa', 'total')+
  theme_classic()+
  scale_fill_tableau(palette='Red-Blue-Brown')+
  labs(y='Плотность SSR \nлокусов (Loci/Mbp)', x='Размер мотива', tag = 'Б')

p3 <- ggplot()+
  geom_boxplot(data=data2, aes(x=tool, y=freq, fill=tool))+
  labs(x='', y='Плотность SSR \nлокусов (Loci/Mbp)', tag = 'В')+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_tableau(palette='Classic Blue-Red 6')+
  theme(plot.margin = unit(c(0.5,0.5,0.5,3), "lines"))



layout <- "
  AAAAA#
  AAAAA#
  AAAAA#
  BBBBCC
  BBBBCC
  "

setwd('G:/Мой диск/Bioinform/Dissertation/Графики')

p <- p1 + p2 + p3 + plot_layout(design = layout) & theme(axis.text=element_text(size=13), axis.title=element_text(size=15),
                                                         legend.text = element_text(size = 13), legend.title = element_text(size=15))

ggplot2::ggsave(filename = "7.Repeats.svg", 
                plot = p, 
                device = "svg", 
                width = 12,
                height = 8, 
                units = "in")



