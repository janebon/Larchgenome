# To make triple plot for paper

library(Biostrings)
library(seqinr)
library(stringr)
library(data.table)

library(ggthemes)
library(gridExtra)
library(ggsci)
library(ggplot2)

library(signs)
library(scales)

# Load data ----------------------------

# Count GC3 of genes  (CDS/intronless, ATG-centered)
setwd('C:/Users/Пользователь/Google диск/Bioinform/TSS/data/For_TRANSFAC/Intronless_genes')
a1 <- 'maker-res.intronless_genes.aligned.fa'
a2 <- 'Pabies1.0.intronless_genes.aligned.fa'
a3 <- 'PG29v3.intronless_genes.aligned.fa'
a4 <- 'Pita.2_01.intronless_genes.aligned.fa'
a5 <- 'TAIR10_cds_20110103_representative_gene_model_updated.txt'  # arabid
a6 <- 'Rice_-_GCF_001433935.1_IRGSP-1.0_cds_from_genomic.fna'
f1 <- read.fasta(a1)
f2 <- read.fasta(a2)
f3 <- read.fasta(a3)
f4 <- read.fasta(a4)
f5 <- read.fasta(a5)
f6 <- read.fasta(a6)
r1 <- as.data.frame(names(f1))
r1$len <- getLength(f1)
r1 <- cbind(r1, as.data.frame(sapply(f1, FUN=GC3)))
names(r1)[3] <- 'gc3'
r2 <- as.data.frame(names(f2))
r2$len <- getLength(f2)
r2 <- cbind(r2, as.data.frame(sapply(f2, FUN=GC3)))
names(r2)[3] <- 'gc3'
r3 <- as.data.frame(names(f3))
r3$len <- getLength(f3)
r3 <- cbind(r3, as.data.frame(sapply(f3, FUN=GC3)))
names(r3)[3] <- 'gc3'
r4 <- as.data.frame(names(f4))
r4$len <- getLength(f4)
r4 <- cbind(r4, as.data.frame(sapply(f4, FUN=GC3)))
names(r4)[3] <- 'gc3'
r5 <- as.data.frame(names(f5))
r5$len <- getLength(f5)
r5 <- cbind(r5, as.data.frame(sapply(f5, FUN=GC3)))
names(r5)[3] <- 'gc3'
r6 <- as.data.frame(names(f6))
r6$len <- getLength(f6)
r6 <- cbind(r6, as.data.frame(sapply(f6, FUN=GC3)))
names(r6)[3] <- 'gc3'

# read GC-gradient values
setwd('C:/Users/Пользователь/Google диск/Bioinform/TSS/data/For_TRANSFAC/Intronless_genes/Computed_data')
c1 <- fread("maker-res.intronless_genes.aligned.fa.FREQS.GC3.TXT")
c2 <- fread("Pabies1.0.intronless_genes.aligned.fa.FREQS.GC3.TXT")
c3 <- fread("PG29v3.intronless_genes.aligned.fa.FREQS.GC3.TXT")
c4 <- fread("Pita.2_01.intronless_genes.aligned.fa.FREQS.GC3.TXT")
c5 <- fread("TAIR10_cds_20110103_representative_gene_model_updated.txt.FREQS.GC3.TXT")  # arabid
c6 <- fread("rice.FREQS.GC3.TXT")

l1 <- lm(GC3_count ~ pos, data = c1, model=TRUE)         # Compute linear regression for GC-gradient
l2 <- lm(GC3_count ~ pos, data = c2, model=TRUE) 
l3 <- lm(GC3_count ~ pos, data = c3, model=TRUE)
l4 <- lm(GC3_count ~ pos, data = c4, model=TRUE)
l5 <- lm(GC3_count ~ pos, data = c5, model=TRUE)
l6 <- lm(GC3_count ~ pos, data = c6, model=TRUE)

coef <- as.data.frame(rbind(l1$coefficients, l2$coefficients, l3$coefficients, 
      l4$coefficients, l5$coefficients, l6$coefficients))    # pull all coefficents together

fu <- cbind(c1$GC3_count, c2$GC3_count, c3$GC3_count, 
            c4$GC3_count, c5$GC3_count, c6$GC3_count)        # pull together GC# counts

coef$gc <- apply(fu, 2, mean)*100                            # calculate mean GC#

coef$name <- c('L.sibirica', 'P.abies', 'P.glauca', 
               'P.taeda', 'A.thaliana', 'O.sativa')                # add names

#      read CG-skew values
setwd('C:/Users/Пользователь/Google диск/Bioinform/TSS/data/For_TRANSFAC/TSS-centered_promoters_new_1000')
s1 <- fread("Larix.TSS-centered.new.1000.fa.FREQS.SKEW.new.TXT")
s2 <- fread("Pabies.TSS-centered.new.1000.fa.FREQS.SKEW.new.TXT")
s3 <- fread("Pglauca.TSS-centered.new.1000.fa.FREQS.SKEW.new.TXT")
s4 <- fread("Pita.TSS-centered.new.1000.fa.FREQS.SKEW.new.TXT")
s5 <- fread("promoters_arabidopsis.fasta.FREQS.SKEW.new.TXT")  # arabidopsis from TT

# Make Plot ------------------------


pal <- c( 'L. sibirica' = "#E64B35FF", 'P. abies' = "#4DBBD5FF", 'P. glauca' = "#00A087FF", 'P. taeda' = "#3C5488FF", 
          'A. thaliana' = "grey60", 'O. sativa' = "#F39B7FFF")
lin <- c( 'L. sibirica'= 1, 'P. abies'=6, 'P. glauca'=5, 'P. taeda'=4, 'A. thaliana'=2, 'O. sativa' = 3)
shps <- c('L. sibirica'= 18, 'P. abies'= 15, 'P. glauca'= 19, 'P. taeda'= 17, 'A. thaliana'= 3, 'O. sativa' = 4)

unicode_minus = function(x) sub('^-', '\U2212', format(x))


p1 <- ggplot()+                      # GC-gradient plot
geom_point(data=c1, aes(x=pos, y=GC3_count, color = 'L. sibirica', shape="L. sibirica"), size=3, alpha=0.5)+
geom_point(data=c2, aes(x=pos, y=GC3_count, color = 'P. abies', shape="P. abies"), size=2, alpha=0.5)+
geom_point(data=c3, aes(x=pos, y=GC3_count, color = 'P. glauca', shape="P. glauca"), size=2, alpha=0.5)+
geom_point(data=c4, aes(x=pos, y=GC3_count, color = 'P. taeda', shape="P. taeda"), size=2, alpha=0.5)+
geom_point(data=c5, aes(x=pos, y=GC3_count, color = 'A. thaliana', shape="A. thaliana"), size=2)+
geom_point(data=c6, aes(x=pos, y=GC3_count, color = 'O. sativa', shape="O. sativa"), size=2)+
  scale_x_continuous(limits = c(1, 1000))+
  scale_y_continuous(limits = c(0.075, 0.45))+
labs(tag = 'a', x = 'Distance from ATG, bp', y='GC3')+
scale_color_manual(name = 'Species', values=pal)+
scale_shape_manual(name='Species', values=shps)+
theme_few()+
theme(axis.text=element_text(size=15), 
      axis.title=element_text(size=15),
      plot.title = element_text(size=18))+
theme(legend.position = c(.97, .97),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    legend.text = element_text(face='italic'))

p2 <- ggplot()+                      # GC3 plot
geom_freqpoly(data=r1, aes(x=gc3, color = 'L. sibirica', linetype='L. sibirica'), size = 1)+
geom_freqpoly(data=r2, aes(x=gc3, color = 'P. abies', linetype='P. abies'), size = 1)+
geom_freqpoly(data=r3, aes(x=gc3, color = 'P. glauca', linetype='P. glauca'), size = 1)+
geom_freqpoly(data=r4, aes(x=gc3, color = 'P. taeda', linetype='P. taeda'), size = 1)+
geom_freqpoly(data=r5, aes(x=gc3, color = 'A. thaliana', linetype='A. thaliana'), size = 1)+
geom_freqpoly(data=r6, aes(x=gc3, color = 'O. sativa', linetype='O. sativa'), size = 1.5)+
labs(tag = 'c', x='GC3', y = 'Number of CDS')+
scale_linetype_manual(name = 'Species', values=lin)+
scale_color_manual(name = 'Species', values=pal)+
theme_few()+
theme(axis.text=element_text(size=15), 
      axis.title=element_text(size=15),
      plot.title = element_text(size=18))+
  theme(legend.position = c(.97, .97),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    legend.text = element_text(face='italic'))


p4 <- ggplot(coef, aes(gc, pos))+           # Slope plot
  geom_point(size=5, color=c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "gray60", "#F39B7FFF"))+
  geom_text(aes(label=name), hjust=0.5, vjust=-2)+
  scale_y_continuous(limits=c(-0.000225, -0.00010), labels = unicode_minus)+
  xlim(c(16,24))+
  labs(tag = 'b', x='GC3', y='Slope')+
  theme_few()


pal <- c( 'L. sibirica' = "#E64B35FF", 'P. abies' = "#4DBBD5FF", 'P. glauca' = "#00A087FF", 'P. taeda' = "#3C5488FF", 'A. thaliana' = "grey60")
lin <- c( 'L. sibirica'= 1, 'P. abies'=6, 'P. glauca'=5, 'P. taeda'=4, 'A. thaliana'=2)
shps <- c('L. sibirica'= 18, 'P. abies'= 15, 'P. glauca'= 19, 'P. taeda'= 17, 'A. thaliana'= 3)
p3 <- ggplot()+                     # Skew plot
  geom_line(data=s1, aes(pos, skew, color = 'L. sibirica', linetype='L. sibirica'), size = 1)+
  geom_line(data=s2, aes(pos, skew, color = 'P. abies', linetype='P. abies'), size = 1)+
  geom_line(data=s3, aes(pos, skew, color = 'P. glauca', linetype='P. glauca'), size = 1)+
  geom_line(data=s4, aes(pos, skew, color = 'P. taeda', linetype='P. taeda'), size = 1)+
  geom_line(data=s5, aes(pos, skew, color = 'A. thaliana', linetype='A. thaliana'), size = 1)+
  scale_x_continuous(limits = c(-500, 500), breaks = scales::pretty_breaks(n = 10), labels = unicode_minus)+
  scale_y_continuous(labels = unicode_minus)+
  geom_vline(xintercept = 0)+ 
  geom_hline(yintercept = 0)+
  scale_color_manual(name = 'Species', values = pal)+
  scale_linetype_manual(name = 'Species', values = lin)+
  labs(tag='d', y='CG-Skew', x='Distance from TSS, bp')+
  theme_few()+
  theme(axis.text=element_text(size=15), 
        axis.title=element_text(size=15),
        plot.title = element_text(size=18))+
  theme(legend.position = c(.97, .97),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.text = element_text(face='italic'))

  setwd('C:/Users/Пользователь/Google диск/Bioinform/TSS/Conifer regulatory regions/Publication')

  png(paste("Fig.6.GC3-Gradient-Skew.png", sep=""), width = 14, height = 9, units = 'in', res = 800)
  
  grid.arrange(p1, p4, p2, p3, layout_matrix = rbind(c(1, 1, 1, 2), c(3, 3, 4, 4)))
  
  dev.off()
