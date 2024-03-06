library(ggsci)
library(gridExtra)
library(ggthemes)
library(data.table)
library(ggplot2)


# Load files
setwd('G:/Мой диск/Bioinform/Dissertation/Доп.мат/Plots/energy')     # here new files are used, from TSS_centered_new_1000, generated on 13.01.2021
f1 <- fread('Larix_stb.txt', header=F, sep='\t', fill=T) 
f2 <- fread('Pabies_stb.txt', header=F, sep='\t', fill=T)
f3 <- fread('Pglauca_stb.txt', header=F, sep='\t', fill=T)
f4 <- read.table('Pita_stb.txt', header=F, fill=T)

# Compute average across all sequences in all windows
avs1 <- as.data.frame(sapply(f1[,-1], 2, FUN=mean, na.rm=T))
names(avs1) <- 'aver_E'
avs1$window <- seq(-999, nrow(avs1)-1000)
avs2 <- as.data.frame(sapply(f2[,-1], 2, FUN=mean, na.rm=T))
names(avs2) <- 'aver_E'
avs2$window <- seq(-999, nrow(avs2)-1000)
avs3 <- as.data.frame(sapply(f3[,-1], 2, FUN=mean, na.rm=T))
names(avs3) <- 'aver_E'
avs3$window <- seq(-999, nrow(avs3)-1000)
avs4 <- as.data.frame(sapply(f4[,-1], 2, FUN=mean, na.rm=T))
names(avs4) <- 'aver_E'
avs4$window <- seq(-999, nrow(avs4)-1000)



pal <- c( 'L. sibirica' = "#E64B35FF", 'P. abies' = "#4DBBD5FF", 'P. glauca' = "#00A087FF", 'P. taeda' = "#3C5488FF")
lin <- c( 'L. sibirica'= 1, 'P. abies'=1, 'P. glauca'=1, 'P. taeda'=1)

unicode_minus = function(x) sub('^-', '\U2212', format(x))

p <- ggplot()+
  geom_line(data=avs1, aes(x=window , y=aver_E, color = 'L. sibirica', linetype='L. sibirica'), size = .7)+
  geom_line(data=avs2, aes(x=window , y=aver_E, color = 'P. abies', linetype='P. abies'), size = .7)+
  geom_line(data=avs3, aes(x=window , y=aver_E, color = 'P. glauca', linetype='P. glauca'), size = .7)+
  geom_line(data=avs4, aes(x=window , y=aver_E, color = 'P. taeda', linetype='P. taeda'), size = .7)+
  geom_vline(xintercept=0)+
  
  scale_color_manual(name='Вид', values=pal)+
  scale_linetype_manual(name = 'Вид', values=lin)+
  scale_x_continuous(lim=c(-900,900), breaks = scales::pretty_breaks(n = 20), labels = unicode_minus) +
  scale_y_continuous(labels = unicode_minus)+
  
  labs(y='Средняя свободная энергия, ккал/моль', x='Расстояние от сайта начала транскрипции, п.н.')+
  
  theme_few()+
  theme(legend.position = c(.85, .8),
        legend.text = element_text(face='italic'))


setwd('G:/Мой диск/Bioinform/Dissertation/Графики')

ggplot2::ggsave(filename = "19.Energy.svg", 
                plot = p, 
                device = "svg", 
                width = 8,
                height = 4, 
                units = "in")

    
