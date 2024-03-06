

setwd('G:/Мой диск/Bioinform/Lar Annot/Repeats/Repeats_to_genes')

library(ggplot2)
library(ggsci)
library(ggthemes)
library(patchwork)

library(dplyr)

getAttributeField <- function (x, field, attrsep = ";") {
  s = strsplit(x$attribute, split = attrsep, fixed = TRUE)
  sapply(s, function(atts) {
    a = strsplit(atts, split = "=", fixed = TRUE)
    m = match(field, sapply(a, "[", 1))
    if (!is.na(m)) {
      rv = a[[m]][2]
    }
    else {
      rv = as.character(NA)
    }
    return(rv)
  })
}


# Count ----------------------------------------------------------------

dt <- read.table('intersect_count.bed2')
names(dt) <- c('chr', 'start', 'stop', 'name', 'score', 'strand', 'count')

rep <- subset(dt, count > 0)  # CDS that intersect with repeats

length(unique(rep$name))      # Number of genes whose CDS intersect with repeats



# wawb  ----------------------------------------------------------------

dt <- read.table('intersect_wawb.bed2')
names(dt) <- c('chrG', 'startG', 'stopG', 'nameG', 'sroceG', 'strandG', 'chrR', 'startR', 'stopR', 'nameR')

tbl <- as.data.frame(table(dt$nameR))
write.table(tbl, 'summary_genes_to_repeats3', quote = FALSE, row.names = FALSE, sep = '\t')


# after manual fix of the table
tbl <- read.table('summary_genes_to_repeats3', header=TRUE)


ggplot(tbl, aes(x = Name, y = Freq, fill = Fam))+
  geom_bar(stat = 'identity', color = 'black', size = 1)+
  facet_wrap(~Fam, scales = 'free_x')+
  labs(y = 'Number of repeats')+
  theme_bw()+
  scale_fill_gdocs()+
  theme(text = element_text(size = 30),
        axis.text.x = element_text(size = 17, angle = 20, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = 'none',
        strip.background = element_rect(fill="grey90"))






## Merge with GFF  -------------------------------------------------------------

gff <- read.table("G:/Мой диск/Bioinform/Lar Annot/maker-res.wi.b2g.incmpl.renamed.NCBI.standardized.bact.gff", sep='\t', quote = "", comment.char = '#')
names(gff) <- c('chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attribute')
mrna <- subset(gff, type == 'mRNA')
mrna$name <- getAttributeField(mrna, 'ID')


dt <- read.table('intersect_count.bed2')
names(dt) <- c('chr', 'start', 'stop', 'name', 'score', 'strand', 'count')
rep <- subset(dt, count > 0)  # CDS that intersect with repeats

length(unique(subset(dt, count > 0)[,'name']))      # Number of genes whose CDS intersect with repeats


dt$name <- gsub(':cds', '', dt$name)
dt2 <- dt %>% group_by(name) %>% top_n(1, count) %>% ungroup()                  # select top count by name (to not drop non-zero counts accidentally)
dt2 <- dt2[!duplicated(dt2[,'name']),]                                          # get one name per gene
dt2 <- dt2[, c('chr', 'score', 'strand', 'name', 'count')]
res <- merge(mrna, dt2)                                                         # merge mRNAs from GFF and repeat counts


# Add Repeat overlap info to GFF

res$att <- ifelse(res$count > 0, paste(res$attribute, ';Note=Overlaps with repeat'), res$attribute)
res$attribute <- NULL ; res$count <- NULL ; res$name <- NULL
names(res)[9] <- 'attribute'

gff2 <- rbind(subset(gff, type!='mRNA'), res)

write.table(gff2, 'G:/Мой диск/Bioinform/Lar Annot/maker-res.wi.b2g.incmpl.renamed.NCBI.standardized.bact.rep.gff', quote = FALSE, row.names = FALSE, sep = '\t')






# Summarise functional categories of genes that overlap with repeats -----------

gff <- read.table("G:/Мой диск/Bioinform/Lar Annot/maker-res.wi.b2g.incmpl.renamed.NCBI.standardized.bact.gff", sep='\t', quote = "", comment.char = '#')
names(gff) <- c('chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attribute')
mrna <- subset(gff, type == 'mRNA')
mrna$name <- getAttributeField(mrna, 'ID')

dt <- read.table('intersect_count.bed2')
names(dt) <- c('chr', 'start', 'stop', 'name', 'score', 'strand', 'count')
dt$name <- gsub(':cds', '', dt$name)
dt2 <- dt %>% group_by(name) %>% top_n(1, count) %>% ungroup()                  # select top count by name (to not drop non-zero counts accidentally)
dt2 <- dt2[!duplicated(dt2[,'name']),]                                          # get one name per gene
dt2 <- dt2[, c('chr', 'score', 'strand', 'name', 'count')]
res <- merge(mrna, dt2)  

rm(dt, dt2, mrna)

repeats <- subset(res, count > 0)
nonsense <- repeats[!grepl("Ontology", repeats$attribute),]
sense <- repeats[grepl("Ontology", repeats$attribute),]
sense$descr <- getAttributeField(sense, 'description')
fu <- as.data.frame(table(sense$descr))



su <- data.frame(fun = c("receptor-like protein kinase", "glutathione S-transferase", 
                   "ABC transporter", "histone ", "SKP1",
                   "glucosyltransferase", "LRR", "cytochrome", "lectin-domain",
                   "pectinesterase", "expansin", "laccase",  "benzoyltransferase", 
                   "abscisic acid",  "peroxidase", "transcription factor", "amino acids", "endoglucanase", 'auxin',
                   "heat shock","aquaporin", "resistance","synthase",
                   "protease", "hydrolase","esterase", "zinc finger","mitochondrial", 
                   'reductase', 'dehydrogenase', "hydroxylase", "lipoxygenase"),
          
                 count = c(sum(fu[grepl("receptor-like protein kinase", fu$Var1),][,2]),
                       sum(fu[grepl("glutathione S-transferase", fu$Var1),][,2]),
                       sum(fu[grepl("ABC transporter", fu$Var1),][,2]),
                       sum(fu[grepl("histone ", fu$Var1),][,2]),
                       sum(fu[grepl("SKP1", fu$Var1),][,2]),
                       sum(fu[grepl("glucosyltransferase", fu$Var1),][,2]),
                       sum(fu[grepl("LRR|leucine-rich", fu$Var1),][,2]),
                       sum(fu[grepl("cytochrome", fu$Var1),][,2]),
                       sum(fu[grepl("lectin-domain", fu$Var1),][,2]) ,
                       sum(fu[grepl("pectinesterase", fu$Var1),][,2]),
                       sum(fu[grepl("expansin", fu$Var1),][,2]),
                       sum(fu[grepl("laccase", fu$Var1),][,2]),
                       sum(fu[grepl("benzoyltransferase", fu$Var1),][,2]),
                       sum(fu[grepl("abscisic", fu$Var1),][,2]),
                       sum(fu[grepl("peroxidase", fu$Var1),][,2]),
                       sum(fu[grepl("transcription factor", fu$Var1),][,2]),
                       sum(fu[grepl("amino acid|amino-", fu$Var1),][,2]),
                       sum(fu[grepl("endoglucanase", fu$Var1),][,2]),
                       sum(fu[grepl("auxin", fu$Var1),][,2]),
                       sum(fu[grepl("heat shock", fu$Var1),][,2]),
                       sum(fu[grepl("aquaporin", fu$Var1),][,2]),
                       sum(fu[grepl("resistance", fu$Var1),][,2]),
                       sum(fu[grepl("synthase", fu$Var1),][,2]),
                       sum(fu[grepl("protease", fu$Var1),][,2]),
                       sum(fu[grepl("hydrolase", fu$Var1),][,2]),
                       sum(fu[grepl("esterase", fu$Var1),][,2]),
                       sum(fu[grepl("zinc finger", fu$Var1),][,2]),
                       sum(fu[grepl("mitochondrial", fu$Var1),][,2]) ,
                       sum(fu[grepl("reductase", fu$Var1),][,2]) ,
                       sum(fu[grepl("dehydrogenase", fu$Var1),][,2]),
                       sum(fu[grepl("hydroxylase", fu$Var1),][,2]),
                       sum(fu[grepl("lipoxygenase", fu$Var1),][,2])))


write.table(su, 'summary_fun_genes3', quote = FALSE, row.names = FALSE, sep = '\t')

barplot(su$count)


ggplot(su, aes(x = reorder(fun, -count), y=count))+
  geom_bar(stat = 'identity', fill='cornflowerblue', color='black')+
  labs(y = 'Number of genes', x = '')+
  theme_bw()+
  theme(axis.text.x = element_text(angle=35, hjust=1),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))


#####--------------------------------------------------------------------------





# Combine plots

tbl <- read.table('summary_genes_to_repeats3', header=TRUE)
su <- read.table('summary_fun_genes3', header=TRUE, sep = '\t')

cols <- c("#CF0811", "#497BA9", "#BA7C99", "#ADB6B6")
  colorRampPalette(pal_lancet()(8))(17)



p1 <- ggplot(tbl, aes(x = Name, y = Freq, fill = Fam))+
  geom_bar(stat = 'identity', color = 'black', size = 1)+
  facet_wrap(~Fam, scales = 'free_x', nrow = 1, ncol  = 4)+
  labs(y = 'Number of repeats')+
  theme_few()+
  scale_fill_manual(values = cols)+
  theme(text = element_text(size = 15),
        axis.text.x = element_text(size = 13, angle = 20, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = 'none')

p2 <- ggplot(su, aes(x = reorder(fun, -count), y=count))+
  geom_bar(stat = 'identity', fill='#b5c8e2', color='black')+
  geom_text(aes(label = fun), angle = 20, hjust = .05, vjust = -1) +
  labs(y = 'Number of genes', x = '')+
  scale_x_discrete(expand = expansion(mult = c(0.02, 0.1)))+
  scale_y_continuous(limits = c(0, 400))+
  theme_few()+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.text.y =  element_text(size = 15))



png("Repeats_to_genes.png" ,  width = 13, height = 10, units = 'in', res = 600)

p1 / p2 + plot_annotation(tag_levels = 'A') + plot_layout(nrow = 2,
                                                          byrow = TRUE,heights = c(1,1.5))


dev.off()














