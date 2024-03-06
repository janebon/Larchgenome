setwd('G:/Мой диск/Bioinform/Lar Annot/B2G/parsing_GO/')

library(ggplot2)
library(patchwork)
library(ggthemes)
library(dplyr)

library(ggsignif)

## Read data --------------
CW <- read.table('G:/Мой диск/Bioinform/Lar Annot/B2G/parsing_GO/For_plots_CW-CD-PCD/CW.txt', header = TRUE, sep = '\t')
CW <- CW %>%  mutate(variable = factor(variable, 
                                       levels = c("PM", "PT", "PL", "PG", "PA", "LS",
                                                  "BP", "PR", "VV", "QR", "FS")))

CD <- read.table('G:/Мой диск/Bioinform/Lar Annot/B2G/parsing_GO/For_plots_CW-CD-PCD/CD.txt', header = TRUE, sep = '\t')
CD <- CD %>%  mutate(variable = factor(variable, 
                                       levels = c("PM", "PT", "PL", "PG", "PA", "LS",
                                                  "BP", "PR", "VV", "QR", "FS")))

HM <- read.table('G:/Мой диск/Bioinform/Lar Annot/B2G/parsing_GO/For_plots_CW-CD-PCD/HM.txt', header = TRUE, sep = '\t')
HM <- HM %>%  mutate(variable = factor(variable, 
                                       levels = c("PM", "PT", "PL", "PG", "PA", "LS",
                                                  "BP", "PR", "VV", "QR", "FS")))

## PALLETES -----------

CW_pal <- c(LS = 'gold',
            PM = 'ghostwhite', 
            PT = 'ghostwhite',
            PL = 'ghostwhite',
            PG = 'ghostwhite',
            PA = 'ghostwhite',
            BP = 'gray30',
            PR = 'gray30',
            VV = 'gray30',
            QR = 'gray30',
            FS = 'gray30' ,
            conifer='ghostwhite', flowering = 'gray30')

CD_pal <- c(LS = 'gold',
            PM = 'ghostwhite', 
            PT = 'ghostwhite',
            PL = 'ghostwhite',
            PG = 'ghostwhite',
            PA = 'ghostwhite',
            BP = 'gray30',
            PR = 'gray30',
            VV = 'gray30',
            QR = 'gray30',
            FS = 'gray30' ,
            conifer='ghostwhite', flowering = 'gray30')


HM_pal <- c(LS = 'gold',
            PM = 'ghostwhite', 
            PT = 'ghostwhite',
            PL = 'ghostwhite',
            PG = 'ghostwhite',
            PA = 'ghostwhite',
            BP = 'gray30',
            PR = 'gray30',
            VV = 'gray30',
            QR = 'gray30',
            FS = 'gray30' ,
            conifer='ghostwhite', flowering = 'gray30')


## Plots -----------------
setwd('G:/Мой диск/Bioinform/Lar Annot/Тексты/Figures/')

## Cell wall ----


p1 <- ggplot(data=CW, aes(x = GO.Names, y = value, fill = variable, width=.5)) +
  geom_bar(stat = 'identity', position = position_dodge(width = 1.2), alpha = 1, color= 'black')+
  geom_text(aes(label = variable), vjust =-0.5, hjust = 0.5, position = position_dodge(1.2), size = 3, angle = 0, color = 'black')+
  facet_wrap(.~ GO.Names, scales="free", strip.position = "bottom", labeller = label_wrap_gen(width=30))+
  labs(x='', y='Percentage of genes', fill='Species')+
  scale_y_continuous(expand = expansion(add = c(0, 0.01)))+
  scale_fill_manual(values = CW_pal)+
  theme_classic() +
  theme(legend.position = 'none', axis.text.x = element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank())

p2 <-  fu 


ggplot() +
  geom_boxplot(data=CW, aes(x=GO.Names, y=value, fill = cat), alpha = 0.8, color = 'black')+
  
  geom_signif(data = CW, aes(x=cat, annotations = c("flowering", "conifer")), manual = TRUE,
    comparisons = list(c('flowering', 'conifer')), map_signif_level = TRUE, textsize = 6) +
  
  facet_wrap(.~ GO.Names, scales="free", strip.position = "bottom", labeller = label_wrap_gen(width=30))+
  labs(x='', y='Percentage of genes', fill='Species')+
  scale_fill_manual(values = CW_pal)+
  theme_classic() +
  theme(legend.position = 'none', axis.text.x = element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank()) 



p <- p1 / p2 & theme(axis.text.y = element_text(size=13), axis.title.y = element_text(size = 20, face = 'bold'), 
                strip.text = element_text(size = 9, face = 'bold'))



ggplot2::ggsave(filename = "Figure-6.GO_Cell_wall.svg", 
                plot = p, 
                device = "svg", 
                width = 10,
                height = 10, 
                units = "in")

# to edit and save as pdf in InkScape



## Programmed Cell death ----

p1 <- ggplot(data=CD, aes(x = GO.Names, y = value, fill = variable, width=.5)) +
  geom_bar(stat = 'identity', position = position_dodge(width = 1.2), alpha = 1, color= 'black')+
  
  geom_text(aes(label = variable), vjust =-0.5, hjust = 0.5, position = position_dodge(1.2), size = 3, angle = 0, color = 'black')+
  
  facet_wrap(.~ GO.Names, scales="free", strip.position = "bottom", ncol = 4, labeller = label_wrap_gen(width=30))+
  labs(x='', y='Percentage of genes', fill='Species')+
  scale_y_continuous(expand = expansion(add = c(0, 0.01)))+
  scale_fill_manual(values = CD_pal)+
  theme_classic() +
  theme(legend.position = 'none', axis.text.x = element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank())

p2 <-  ggplot() +
  geom_boxplot(data=CD, aes(x=GO.Names, y=value, fill = cat), alpha = 0.8, color = 'black')+
  
  facet_wrap(.~ GO.Names, scales="free", strip.position = "bottom", ncol = 4, labeller = label_wrap_gen(width=30))+
  labs(x='', y='Percentage of genes', fill='Species')+
  scale_fill_manual(values = CD_pal)+
  theme_classic() +
  theme(legend.position = 'none', axis.text.x = element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(),) 


p <- p1 / p2 & theme(axis.text.y = element_text(size=13), axis.title.y = element_text(size = 20, face = 'bold'), 
                     strip.text = element_text(size = 9, face = 'bold'))


ggplot2::ggsave(filename = "Figure-7.GO_PCD.svg", 
                plot = p, 
                device = "svg", 
                width = 10,
                height = 4, 
                units = "in")

# to edit and save as pdf in InkScape






## Hormones -----

p1 <- ggplot(data=HM, aes(x = GO.Names, y = value, fill = variable, width=.5)) +
  geom_bar(stat = 'identity', position = position_dodge(width = 1.2), alpha = 1, color= 'black')+
  
  geom_text(aes(label = variable), vjust =-0.5, hjust = 0.5, position = position_dodge(1.2), size = 3, angle = 0, color = 'black')+
  
  facet_wrap(.~ GO.Names, scales="free", strip.position = "bottom", ncol = 4, labeller = label_wrap_gen(width=30))+
  labs(x='', y='Percentage of genes', fill='Species')+
  scale_y_continuous(expand = expansion(add = c(0, 0.02)))+
  scale_fill_manual(values = HM_pal)+
  theme_classic() +
  theme(legend.position = 'none', axis.text.x = element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank())

p2 <-  ggplot() +
  geom_boxplot(data=HM, aes(x=GO.Names, y=value, fill = cat), alpha = 0.8, color = 'black')+
  
  facet_wrap(.~ GO.Names, scales="free", strip.position = "bottom", ncol = 4, labeller = label_wrap_gen(width=30))+
  labs(x='', y='Percentage of genes', fill='Species')+
  scale_fill_manual(values = HM_pal)+
  theme_classic() +
  theme(legend.position = 'none', axis.text.x = element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank())


p <- p1 / p2 & theme(axis.text.y = element_text(size=13), axis.title.y = element_text(size = 20, face = 'bold'), 
                     strip.text = element_text(size = 9, face = 'bold'))

ggplot2::ggsave(filename = "Figure-8.GO_Hormones.svg", 
                plot = p, 
                device = "svg", 
                width = 10,
                height = 7, 
                units = "in")

# to edit and save as pdf in InkScape


