
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

cols = c('Lar' = "#3366cc", 'Pila' = "#dc3912", 'Pita' = "#ff9900", 'Pab' = "#109618", 'Pgl' = "#990099", 'Psme' = "#66aa00", 
         'Fasy' = "#0099c6", 'Quro' = "#dd4477", 'Potr' = "#b82e2e", 'Vivi' = "#316395")

lns = c('Lar' = 1, 'Pila' = 2, 'Pita' = 2, 'Pab' = 2, 'Pgl' = 2, 'Psme' = 2, 
         'Fasy' = 3, 'Quro' = 3, 'Potr' = 3, 'Vivi' = 3)

  
p2 <- ggplot() + 
  stat_ecdf(data=Lar, aes(-qcovhsp, color = 'Lar', linetype = 'Lar'), geom = "line", size=1)+
  
  stat_ecdf(data=Pila, aes(-qcovhsp, color = 'Pila', linetype = 'Pila'), geom = "line", size=1)+
  stat_ecdf(data=Pita, aes(-qcovhsp, color = 'Pita', linetype = 'Pita'), geom = "line", size=1)+
  stat_ecdf(data=Pab, aes(-qcovhsp, color = 'Pab', linetype = 'Pab'), geom = "line", size=1)+
  stat_ecdf(data=Pgl, aes(-qcovhsp, color = 'Pgl', linetype = 'Pgl'), geom = "line", size=1)+
  stat_ecdf(data=Psme, aes(-qcovhsp, color = 'Psme', linetype = 'Psme'), geom = "line", size=1)+
  
  stat_ecdf(data=Fasy, aes(-qcovhsp, color = 'Fasy', linetype = 'Fasy'), geom = "line", size=1)+
  stat_ecdf(data=Quro, aes(-qcovhsp, color = 'Quro', linetype = 'Quro'), geom = "line", size=1)+
  stat_ecdf(data=Potr, aes(-qcovhsp, color = 'Potr', linetype = 'Potr'), geom = "line", size=1)+
  stat_ecdf(data=Vivi, aes(-qcovhsp, color = 'Vivi', linetype = 'Vivi'), geom = "line", size=1)+
  
  ylim(0.0001,0.999)+
  scale_x_continuous(breaks = seq(-100, -10, by=10), labels= c('>100','>90','>80','>70','>60','>50','>40',' >30','>20',' >10'))+
  labs(x = 'Query coverage, %', y = 'Proportion of genes')+
  scale_linetype_manual(values = lns, name = 'Species')+
  scale_color_manual(values = cols, name = 'Species')+
  theme_bw()


# Plot NUMBER of genes vs qcovhsp%, cumulative

p1 <- ggplot()+
  geom_point(data=dt.Lar, aes(x=-cutoff, y=Ngenes, color = 'Lar'), pch=1, size=3)+
  geom_line(data=dt.Lar, aes(x=-cutoff, y=Ngenes, color = 'Lar'), size=1)+
  
    geom_point(data=dt.Pila, aes(x=-cutoff, y=Ngenes, color = 'Pila'), pch=1, size=3)+
  geom_line(data=dt.Pila, aes(x=-cutoff, y=Ngenes, color = 'Pila'), size=1, linetype = 2)+
    geom_point(data=dt.Pita, aes(x=-cutoff, y=Ngenes, color = 'Pita'), pch=1, size=3)+
  geom_line(data=dt.Pita, aes(x=-cutoff, y=Ngenes, color = 'Pita'), size=1, linetype = 2)+
    geom_point(data=dt.Pab, aes(x=-cutoff, y=Ngenes, color = 'Pab'), pch=1, size=3)+
  geom_line(data=dt.Pab, aes(x=-cutoff, y=Ngenes, color = 'Pab'), size=1, linetype = 2)+
    geom_point(data=dt.Pgl, aes(x=-cutoff, y=Ngenes, color = 'Pgl'), pch=1, size=3)+
  geom_line(data=dt.Pgl, aes(x=-cutoff, y=Ngenes, color = 'Pgl'), size=1, linetype = 2)+
    geom_point(data=dt.Psme, aes(x=-cutoff, y=Ngenes, color = 'Psme'), pch=1, size=3)+
  geom_line(data=dt.Psme, aes(x=-cutoff, y=Ngenes, color = 'Psme'), size=1, linetype = 2)+
  
  geom_point(data=dt.Fasy, aes(x=-cutoff, y=Ngenes, color = 'Fasy'), pch=1, size=3)+
  geom_line(data=dt.Fasy, aes(x=-cutoff, y=Ngenes, color = 'Fasy'), size=1, linetype = 3)+
    geom_point(data=dt.Quro, aes(x=-cutoff, y=Ngenes, color = 'Quro'), pch=1, size=3)+
  geom_line(data=dt.Quro, aes(x=-cutoff, y=Ngenes, color = 'Quro'), size=1, linetype = 3)+
    geom_point(data=dt.Potr, aes(x=-cutoff, y=Ngenes, color = 'Potr'), pch=1, size=3)+
  geom_line(data=dt.Potr, aes(x=-cutoff, y=Ngenes, color = 'Potr'), size=1, linetype = 3)+
    geom_point(data=dt.Vivi, aes(x=-cutoff, y=Ngenes, color = 'Vivi'), pch=1, size=3)+
  geom_line(data=dt.Vivi, aes(x=-cutoff, y=Ngenes, color = 'Vivi'), size=1, linetype = 3)+
  
  scale_x_continuous(breaks = seq(-100, 0, by=10), labels= c('>100','>90','>80','>70','>60','>50','>40',' >30','>20',' >10','0'))+
  scale_color_manual(values = cols, name = 'Species')+
  theme_bw()+
  labs(x = 'Query coverage, %', y = 'Number of genes')+
  guides(color = FALSE)



svg('Cumulative_plot.svg', width = 14, height = 5)

p1 + p2 + plot_annotation(tag_levels = 'A') & theme(axis.text = element_text(size = 17), axis.title = element_text(size = 17),
                                                    legend.title = element_text(size = 17), legend.text = element_text(size = 17),
                                                    panel.grid.minor = element_blank())

dev.off()



png('Cumulative_plot.png', width = 14, height = 5, res = 600, units = 'in')

p1 + p2 + plot_annotation(tag_levels = 'A') & theme(axis.text = element_text(size = 17), axis.title = element_text(size = 17),
                                                    legend.title = element_text(size = 17), legend.text = element_text(size = 17),
                                                    panel.grid.minor = element_blank())

dev.off()
