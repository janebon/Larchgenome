library(Biostrings)
library(stringr)
library(ggsci)
library(ggthemes)
library(data.table)
library(ggplot2)

#        Use already computed values --------------------------------------

setwd('G:/Мой диск/Bioinform/Dissertation/Доп.мат/Plots/tata')
count_a <- fread('Pabies.new.tata-freq.w20s10.txt'  )
count_b <- fread('Pita.new.tata-freq.w20s10.txt'    )
count_c <- fread('P.glauca.new.tata-freq.w20s10.txt')
count_d <- fread('Larix.new.tata-freq.w20s10.txt'   )
count_z <- fread('TAIR10.new.tata-freq.w20s10.txt'  )
count_a$W <- seq(-1000, abs(100-nrow(count_a))*10-1, by=10)
count_a$W <- count_a$W + 20
count_b$W <- seq(-1000, abs(100-nrow(count_b))*10-1, by=10)
count_b$W <- count_b$W + 20
count_c$W <- seq(-1000, abs(100-nrow(count_c))*10-1, by=10)
count_c$W <- count_c$W + 20
count_d$W <- seq(-1000, abs(100-nrow(count_d))*10-1, by=10)
count_d$W <- count_d$W + 20
count_z$W <- seq(-1000, abs(100-nrow(count_z))*10-1, by=10)
count_z$W <- count_z$W + 20

pal <- c( 'L.sibirica' = "#E64B35FF", 'P.abies' = "#4DBBD5FF", 'P.glauca' = "#00A087FF", 'P.taeda' = "#3C5488FF", 'A.thaliana' = "grey60")
unicode_minus = function(x) sub('^-', '\U2212', format(x))

#      Make a plot with all-in-one motif profile 

p.full <- ggplot()+
  geom_line(data=count_z, aes(x = W, y = signals, color='A.thaliana'), size=1)+
  geom_line(data=count_d, aes(x = W, y = signals, color='P.taeda'), size=1)+
  geom_line(data=count_c, aes(x = W, y = signals, color='P.glauca'), size=1)+
  geom_line(data=count_b, aes(x = W, y = signals, color='P.abies'), size=1)+
  geom_line(data=count_a, aes(x = W, y = signals, color='L.sibirica'), size=1)+
  
  scale_x_continuous(limits = c(-600, 600), breaks = scales::pretty_breaks(n = 10), labels = unicode_minus)+
  scale_color_manual(name = 'Вид', values = pal)+
  geom_vline(xintercept=0)+
  
  labs(x='Расстояние от сайта начала транскрипции, п.н.', y='Число сигналов')+
  
  theme_few()+
  theme(axis.text.x=element_text())+
  theme(legend.position = c(.85, .75),
        legend.text = element_text(face='italic'))

p.zoom <- ggplot()+
  geom_line(data=count_z, aes(x = W, y = signals, color='A.thaliana'), size=1)+
  geom_line(data=count_d, aes(x = W, y = signals, color='P.taeda'), size=1)+
  geom_line(data=count_c, aes(x = W, y = signals, color='P.glauca'), size=1)+
  geom_line(data=count_b, aes(x = W, y = signals, color='P.abies'), size=1)+
  geom_line(data=count_a, aes(x = W, y = signals, color='L.sibirica'), size=1)+
  
  scale_x_continuous(limits = c(-50, 20), breaks = scales::pretty_breaks(), labels = unicode_minus)+
  scale_color_manual(name = 'Вид', values = pal)+
  geom_vline(xintercept=0)+
  
  theme_few()+
  theme(axis.text.x=element_text())+
  theme(legend.position = 'none',
        axis.title=element_blank())

p <- p.full + annotation_custom(grob = ggplotGrob(p.zoom), xmin=-650, xmax=-200, ymin=700, ymax=1300)


setwd('G:/Мой диск/Bioinform/Dissertation/Графики')

ggplot2::ggsave(filename = "18.TATA.svg", 
                plot = p, 
                device = "svg", 
                width = 8,
                height = 4, 
                units = "in")



  
