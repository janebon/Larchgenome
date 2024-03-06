library(ggthemes)
library(gridExtra)
library(ggsci)
library(data.table)
library(ggplot2)
library(scales)

library(patchwork)
library(ggsignif)


# Read count GC3 of genes  (CDS/intronless, ATG-centered) -------------------------
setwd('G:/Мой диск/Bioinform/Dissertation/Доп.мат/Plots/gc3')

r1 <- fread('Larch.gc3gc3.txt')
r2 <- fread('Pabies.gc3gc3.txt')
r3 <- fread('PG.gc3gc3.txt')
r4 <- fread('Pita.gc3gc3.txt')
r5 <- fread('TAIR.gc3gc3.txt')
r6 <- fread('Rice.gc3gc3.txt')

#       Devide genes into Poor and Rich based on quantiles 10-90  ------------------

r1$type <- ifelse(r1$gc3 <= quantile(r1$gc3, 0.1), 'poor', 'middle')
r1$type <- ifelse(r1$gc3 >= quantile(r1$gc3, 0.9), 'rich', r1$type)
r1$type <- factor(r1$type, levels=c("poor", "middle", "rich"), ordered = TRUE)

r2$type <- ifelse(r2$gc3 <= quantile(r2$gc3, 0.1), 'poor', 'middle')
r2$type <- ifelse(r2$gc3 >= quantile(r2$gc3, 0.9), 'rich', r2$type)
r2$type <- factor(r2$type, levels=c("poor", "middle", "rich"), ordered = TRUE)

r3$type <- ifelse(r3$gc3 <= quantile(r3$gc3, 0.1), 'poor', 'middle')
r3$type <- ifelse(r3$gc3 >= quantile(r3$gc3, 0.9), 'rich', r3$type)
r3$type <- factor(r3$type, levels=c("poor", "middle", "rich"), ordered = TRUE)

r4$type <- ifelse(r4$gc3 <= quantile(r4$gc3, 0.1), 'poor', 'middle')
r4$type <- ifelse(r4$gc3 >= quantile(r4$gc3, 0.9), 'rich', r4$type)
r4$type <- factor(r4$type, levels=c("poor", "middle", "rich"), ordered = TRUE)

r5$type <- ifelse(r5$gc3 <= quantile(r5$gc3, 0.1), 'poor', 'middle')
r5$type <- ifelse(r5$gc3 >= quantile(r5$gc3, 0.9), 'rich', r5$type)
r5$type <- factor(r5$type, levels=c("poor", "middle", "rich"), ordered = TRUE)

r6$type <- ifelse(r6$gc3 <= quantile(r6$gc3, 0.1), 'poor', 'middle')
r6$type <- ifelse(r6$gc3 >= quantile(r6$gc3, 0.9), 'rich', r6$type)
r6$type <- factor(r6$type, levels=c("poor", "middle", "rich"), ordered = TRUE)

# Remove middle category
r1 <- subset(r1, type != 'middle'); r1$type <- factor(r1$type)
r2 <- subset(r2, type != 'middle'); r2$type <- factor(r2$type)
r3 <- subset(r3, type != 'middle'); r3$type <- factor(r3$type)
r4 <- subset(r4, type != 'middle'); r4$type <- factor(r4$type)
r5 <- subset(r5, type != 'middle'); r5$type <- factor(r5$type)
r6 <- subset(r6, type != 'middle'); r6$type <- factor(r6$type)

# Test for difference ---------

test1 <- wilcox.test(x=as.matrix(subset(r1, type=='poor')[,2]), y=as.matrix(subset(r1, type=='rich')[,2]), # Siberian larch
                       paired=FALSE, alternative = 'two.sided')$p.value
wilcox.test(x=as.matrix(subset(r2, type=='poor')[,2]), y=as.matrix(subset(r2, type=='rich')[,2]), # Norway spruce
                       paired=FALSE, alternative = 'two.sided')
wilcox.test(x=as.matrix(subset(r3, type=='poor')[,2]), y=as.matrix(subset(r3, type=='rich')[,2]), # White spruce
                       paired=FALSE, alternative = 'two.sided')
wilcox.test(x=as.matrix(subset(r4, type=='poor')[,2]), y=as.matrix(subset(r5, type=='rich')[,2]), # Loblolly pine
                       paired=FALSE, alternative = 'two.sided')

# Plot ------------------------

p1 <- ggplot(r1, aes(x=type, y=len, color=type))+
  geom_point(alpha = 0.3, size=0.5, position = "jitter")+
  geom_boxplot(width=0.5, notch = TRUE, notchwidth = 0.5, outlier.alpha = 0, linewidth = 1)+
  geom_signif(comparisons = list(c("poor", "rich")), map_signif_level = TRUE, textsize = 6, y_position = 7000) +
  scale_color_aaas()+
  scale_y_continuous(labels=comma)+
  labs(title = 'L. sibirica', x=NULL, y='Длина гена, п.н.')+
  theme_few()+
  theme(plot.title=element_text(hjust=0.5, face='italic'), legend.position='none',
        panel.background = element_rect(fill='lightsteelblue1'))

p2 <- ggplot(r2, aes(x=type, y=len, color=type))+
  geom_point(alpha = 0.3, size=0.5, position = "jitter")+
  geom_boxplot(width=0.5, notch = TRUE, notchwidth = 0.5, outlier.alpha = 0, linewidth = 1)+
  geom_signif(comparisons = list(c("poor", "rich")), map_signif_level = TRUE, textsize = 6, y_position = 6500) +
  scale_color_aaas()+
  scale_y_continuous(labels=comma)+
  labs(title = 'P. abies', x=NULL, y='Длина гена, п.н.')+
  theme_few()+
  theme(plot.title=element_text(hjust=0.5, face='italic'), legend.position="none",
        panel.background = element_rect(fill='lightsteelblue1'))

p3 <- ggplot(r3, aes(x=type, y=len, color=type))+
  geom_point(alpha = 0.3, size=0.5, position = "jitter")+
  geom_boxplot(width=0.5, notch = TRUE, notchwidth = 0.5, outlier.alpha = 0, linewidth = 1)+
  scale_color_aaas()+
  scale_y_continuous(labels=comma)+
  labs(title = 'P. glauca', x=NULL, y='Длина гена, п.н.')+
  theme_few()+
  theme(plot.title=element_text(hjust=0.5, face='italic'), legend.position="none",
        panel.background = element_rect(fill='lightsteelblue1'))

p4 <- ggplot(r4, aes(x=type, y=len, color=type))+
  geom_point(alpha = 0.3, size=0.5, position = "jitter")+
  geom_boxplot(width=0.5, notch = TRUE, notchwidth = 0.5, outlier.alpha = 0, linewidth = 1)+
  geom_signif(comparisons = list(c("poor", "rich")), map_signif_level = TRUE, textsize = 6, y_position = 9000) +
  scale_color_aaas()+
  scale_y_continuous(labels=comma)+
  labs(title = 'P. taeda', x=NULL, y='Длина гена, п.н.')+
  theme_few()+
  theme(plot.title=element_text(hjust=0.5, face='italic'), legend.position="none",
        panel.background = element_rect(fill='lightsteelblue1'))

p5 <- ggplot(r5, aes(x=type, y=len, color=type))+
  geom_point(alpha = 0.3, size=0.5, position = "jitter")+
  geom_boxplot(width=0.5, notch = TRUE, notchwidth = 0.5, outlier.alpha = 0, linewidth = 1)+
  geom_signif(comparisons = list(c("poor", "rich")), map_signif_level = TRUE, textsize = 6, y_position = 12000) +
  scale_colour_aaas()+
  scale_y_continuous(labels=comma)+
  labs(title = 'A. thaliana', x=NULL, y='Длина гена, п.н.')+
  theme_few()+
  theme(plot.title=element_text(hjust=0.5, face='italic'), 
        legend.position="none", 
        panel.background = element_rect(fill='gray90'))

p6 <- ggplot(r6, aes(x=type, y=len, color=type))+
  geom_point(alpha = 0.3, size=0.5, position = "jitter")+
  geom_boxplot(width=0.5, notch = TRUE, notchwidth = 0.5, outlier.alpha = 0, linewidth = 1)+
  geom_signif(comparisons = list(c("poor", "rich")), map_signif_level = TRUE, textsize = 6, y_position = 14000) +
  scale_color_aaas()+
  scale_y_continuous(labels=comma)+
  labs(title = 'O. sativa', x=NULL, y='Длина гена, п.н.')+
  theme_few()+
  theme(plot.title=element_text(hjust=0.5, face='italic'),
        legend.position="none", 
        panel.background = element_rect(fill='gray90'))



setwd('G:/Мой диск/Bioinform/Dissertation/Графики')

p <- (p1+p2)/(p3+p4)/(p5+p6) 

ggplot2::ggsave(filename = "22.Boxplot.svg", 
                plot = p, 
                device = "svg", 
                width = 10,
                height = 8, 
                units = "in")








