
setwd('D:/larch/Cumulative plot')

library(ggplot2)
library(ggthemes)
library(patchwork)
library(dplyr)


FILES=list.files(pattern = '*.psl')

for(i in 1:length(FILES)){
  
  f=FILES[i]
  
  dt <- read.table(f, header = FALSE, sep = '\t', quote = '', comment.char = '')
  names(dt) <- c('qseqid', 'stitle', 'qlen', 'slen', 'length', 'qcovhsp', 'evalue', 'pident')
  
  dt2 <- dt %>% group_by(qseqid) %>% 
  slice_min(order_by = evalue, n = 1) %>%
  slice_min(order_by = pident, n = 1) %>%
  slice_min(order_by = qcovhsp, n = 1) %>%
  slice(1) %>%
  ungroup()
  
  dt2 <- as.data.frame(dt2)
  
  nom <- paste(gsub(".arabid.psl","", FILES[i]))
  assign(nom, dt2)
  
  dt3 <- data.frame(c(90,80,70,60,50,40,30,20,10),
                    c(nrow(subset(dt2, qcovhsp > 90)),
                      nrow(subset(dt2, qcovhsp > 80)),
                      nrow(subset(dt2, qcovhsp > 70)),
                      nrow(subset(dt2, qcovhsp > 60)),
                      nrow(subset(dt2, qcovhsp > 50)),
                      nrow(subset(dt2, qcovhsp > 40)),
                      nrow(subset(dt2, qcovhsp > 30)),
                      nrow(subset(dt2, qcovhsp > 20)),
                      nrow(subset(dt2, qcovhsp > 10))))
  names(dt3) <- c('cutoff', 'Ngenes')
  
  mon <- paste('dt.', gsub(".arabid.psl","", FILES[i]), sep='')
  assign(mon, dt3)
  
}

Lar <- Larch.bact ; rm(Larch.bact)
dt.Lar <- dt.Larch.bact ; rm(dt.Larch.bact)


## Plot PROPORTION of genes vs qcovhsp%, cumulative
ggthemes::ggthemes_data$gdocs$colors$value
colorRampPalette(gdocs_pal()(10))(12)

cols = c('L.sibirica' = "#FF3030", 'P.lambertiana' = "#63B8FF", 'P.taeda' = "#FFA500", 'P.abies' = "#008B45", 'P.glauca' = "#FF3E96", 'P.menziesii' = "#A05B8C",
         'A.alba' = '#837C10', 'P.tabuliformis' = '#32CD32',
         'F.sylvatica' = "#3366cc", 'Q.robus' = "#9B30FF", 'P.trichocarpa' = "#b82e2e", 'V.vinifera' = "#316395")

lns = c('L.sibirica' = 1, 
        'P.lambertiana' = 2, 'P.taeda' = 2, 'P.abies' = 2, 'P.glauca' = 2, 'P.menziesii' = 2, 'A.alba' = 2, 'P.tabuliformis' = 2,
         'F.sylvatica' = 3, 'Q.robus' = 3, 'P.trichocarpa' = 3, 'V.vinifera' = 3)



cols = c('L.sibirica' = "violetred1", 
         'P.lambertiana' = "limegreen", 'P.taeda' = "limegreen", 'P.abies' = "limegreen", 'P.glauca' = "limegreen", 'P.menziesii' = "limegreen",
         'A.alba' = 'limegreen', 'P.tabuliformis' = 'limegreen',
         'F.sylvatica' = "steelblue1", 'Q.robus' = "steelblue1", 'P.trichocarpa' = "steelblue1", 'V.vinifera' = "steelblue1")

lns = c('L.sibirica' = 1, 
        'P.lambertiana' = 1, 'P.taeda' = 1, 'P.abies' = 1, 'P.glauca' = 1, 'P.menziesii' = 1, 'A.alba' = 1, 'P.tabuliformis' = 1,
        'F.sylvatica' = 1, 'Q.robus' = 1, 'P.trichocarpa' = 1, 'V.vinifera' = 1)
  
p2 <- ggplot() + 
  stat_ecdf(data=Lar, aes(-qcovhsp, color = 'L.sibirica', linetype = 'L.sibirica'), geom = "line", size=1.3)+
  
  stat_ecdf(data=Pila, aes(-qcovhsp, color = 'P.lambertiana', linetype = 'P.lambertiana'), geom = "line", size=1)+
  stat_ecdf(data=Pita, aes(-qcovhsp, color = 'P.taeda', linetype = 'P.taeda'), geom = "line", size=1)+
  stat_ecdf(data=Pab, aes(-qcovhsp, color = 'P.abies', linetype = 'P.abies'), geom = "line", size=1)+
  stat_ecdf(data=Pgl, aes(-qcovhsp, color = 'P.glauca', linetype = 'P.glauca'), geom = "line", size=1)+
  stat_ecdf(data=Psme, aes(-qcovhsp, color = 'P.menziesii', linetype = 'P.menziesii'), geom = "line", size=1)+
  stat_ecdf(data=Abal, aes(-qcovhsp, color = 'A.alba', linetype = 'A.alba'), geom = "line", size=1)+
  stat_ecdf(data=P.tabuliformis, aes(-qcovhsp, color = 'P.tabuliformis', linetype = 'P.tabuliformis'), geom = "line", size=1)+
  
  stat_ecdf(data=Fasy, aes(-qcovhsp, color = 'F.sylvatica', linetype = 'F.sylvatica'), geom = "line", size=1)+
  stat_ecdf(data=Quro, aes(-qcovhsp, color = 'Q.robus', linetype = 'Q.robus'), geom = "line", size=1)+
  stat_ecdf(data=Potr, aes(-qcovhsp, color = 'P.trichocarpa', linetype = 'P.trichocarpa'), geom = "line", size=1)+
  stat_ecdf(data=Vivi, aes(-qcovhsp, color = 'V.vinifera', linetype = 'V.vinifera'), geom = "line", size=1)+
  
  ylim(0.0001,0.999)+
  scale_x_continuous(breaks = seq(-100, -10, by=20), labels= c('>100','>80','>60','>40','>20'))+
  labs(x = 'Степень покрытия гена, %', y = 'Доля генов')+
  scale_linetype_manual(values = lns, name = 'Вид')+
  scale_color_manual(values = cols, name = 'Вид')+
  theme_bw()+
  theme(legend.text = element_text(face = 'italic'))


# Plot NUMBER of genes vs qcovhsp%, cumulative

p1 <- ggplot()+
  geom_line(data=dt.Lar, aes(x=-cutoff, y=Ngenes, color = 'L.sibirica', linetype = 'L.sibirica'), size=1)+
  geom_line(data=dt.Pila, aes(x=-cutoff, y=Ngenes, color = 'P.lambertiana', linetype = 'P.lambertiana'), size=1)+
  geom_line(data=dt.Pita, aes(x=-cutoff, y=Ngenes, color = 'P.taeda', linetype = 'P.taeda'), size=1)+
  geom_line(data=dt.Pab, aes(x=-cutoff, y=Ngenes, color = 'P.abies', linetype = 'P.abies'), size=1)+
  geom_line(data=dt.Pgl, aes(x=-cutoff, y=Ngenes, color = 'P.glauca', linetype = 'P.glauca'), size=1)+
  geom_line(data=dt.Psme, aes(x=-cutoff, y=Ngenes, color = 'P.menziesii', linetype = 'P.menziesii'), size=1)+
  
  geom_line(data=dt.Abal, aes(x=-cutoff, y=Ngenes, color = 'A.alba', linetype = 'A.alba'), size=1)+
  geom_line(data=dt.P.tabuliformis, aes(x=-cutoff, y=Ngenes, color = 'P.tabuliformis', linetype = 'P.tabuliformis'), size=1)+
  geom_line(data=dt.Fasy, aes(x=-cutoff, y=Ngenes, color = 'F.sylvatica', linetype = 'F.sylvatica'), size=1)+
  geom_line(data=dt.Quro, aes(x=-cutoff, y=Ngenes, color = 'Q.robus', linetype = 'Q.robus'), size=1)+
  geom_line(data=dt.Potr, aes(x=-cutoff, y=Ngenes, color = 'P.trichocarpa', linetype = 'P.trichocarpa'), size=1)+
  geom_line(data=dt.Vivi, aes(x=-cutoff, y=Ngenes, color = 'V.vinifera', linetype = 'V.vinifera'), size=1)+
  
  scale_x_continuous(breaks = seq(-100, 0, by=20), labels= c('>100','>80','>60','>40','>20', '0'))+
  scale_color_manual(values = cols, name = 'Вид')+
  scale_linetype_manual(values = lns, name = 'Вид')+
  theme_bw()+
  labs(x = 'Степень покрытия гена, %', y = 'Число генов')+
  theme(legend.position="none")






p <- p1 + p2 + plot_annotation(tag_levels = 'A') & theme(axis.text = element_text(size = 18,), 
                                                    axis.title = element_text(size = 20),
                                                    legend.text = element_text(size = 18), legend.title = element_text(size = 20), 
                                                    panel.grid.minor = element_blank(), 
                                                    plot.tag = element_text(size = 20, face = 'bold'))

setwd('G:/Мой диск/Bioinform/Dissertation/Графики')

ggplot2::ggsave(filename = "10.Cumulative_plot-2.jpg", 
                plot = p, 
                device = "jpeg", 
                width = 13,
                height = 5, 
                dpi = 600,
                units = "in")

ggplot2::ggsave(filename = "10.Cumulative_plot-2.svg", 
                plot = p, 
                device = "svg", 
                width = 13,
                height = 5, 
                units = "in")






