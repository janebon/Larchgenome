
setwd('G:/Мой диск/Bioinform/Dissertation/mito')
library(ggplot2)
library(ggrepel)
library(svglite)
library(dplyr)


library(ggthemes)
library(ggsci)
library(RColorBrewer)


tab <- read.table('Genes.txt', header=TRUE, sep = '\t')
tab$molecule <- factor(tab$molecule,
                       levels = c('MT797187.1', 'MT797188.1', 'MT797189.1', 'MT797190.1', 'MT797191.1', 'MT797192.1', 'MT797193.1', 'MT797194.1', 'MT797195.1'))


pal <- colorRampPalette(pal_npg()(7))(11)
pal <- colorRampPalette(brewer.pal(9,"Set3"))(11)  # for brewer pal
pal <- colorRampPalette(few_pal()(8))(11)        # for gdocs pal
pal <- colorRampPalette(tableau_color_pal( palette = "Classic Blue-Red 6")(6))(11)


p <- gg

ggplot(tab, aes(x = as.numeric(start_scaff), xend = as.numeric(end_scaff), y = molecule, yend = molecule)) +
  geom_segment(size = 5, col = "grey80") +
  geom_segment(aes(x = ifelse(direction == 1, start_gene, end_gene), xend = ifelse(direction == 1, end_gene, start_gene), color=group), data = tab, arrow = arrow(type = 'closed', length = unit(0.1, "inches")), size = 1) +
  geom_text_repel(data = tab, aes(x = start_gene, y = molecule, label = gene), size = 3, fontface='bold.italic', nudge_y = 0.4) + 
  geom_text(data = tab, aes(x = end_scaff, y = molecule, label = paste0(end_scaff, ' bp')), hjust=-0.2, size = 4, fontface='bold', color = 'grey50', check_overlap = T) + 
  scale_y_discrete(limits = rev(levels(tab$molecule))) +
  expand_limits(x = c(0, 4400000)) +
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank())+
  theme(axis.title = element_blank(), axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size = 15, face = 'bold'))+
  theme(legend.position = c(.98, .42),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.text = element_text(size = 12))+
  scale_colour_manual(values = pal)

ggplot2::ggsave(filename = "Gene_map.svg", 
                plot = p, 
                device = "svg", 
                width = 9.5,
                height = 8, 
                dpi = 600,
                units = "in")
