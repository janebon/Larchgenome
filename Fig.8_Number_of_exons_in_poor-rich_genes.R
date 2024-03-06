library(ggsci)
library(ggthemes)
library(patchwork)
library(ggplot2)


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

exons_merge<- function(gff, gc3){
  exons <- subset(gff, type=='exon')
  exons$name <- getAttributeField(exons, field='Parent')
  exons$name <- gsub("-mRNA-1", "", exons$name)
  exons$name <- gsub(".mrna1", "", exons$name)
  
  poor <- merge(exons, subset(gc3, type=='poor'), by='name')
  rich <- merge(exons, subset(gc3, type=='rich'), by='name')
  
  result <- rbind(poor, rich)
  result <- result[, -4]
  names(result)[12] <- 'type'
  return(result)
}       #  Function to merge gff and table with GC3 and poor-rich types


# Read count GC3 of genes  (CDS/intronless, ATG-centered)  ----------------
  setwd('G:/Мой диск/Bioinform/Dissertation/Доп.мат/Plots/gc3')
r1 <- read.table('Larch.gc3gc3.txt', header = TRUE)
r2 <- read.table('Pabies.gc3gc3.txt', header = TRUE)
r3 <- read.table('PG.gc3gc3.txt', header = TRUE)
r4 <- read.table('Pita.gc3gc3.txt', header = TRUE)
r5 <- read.table('TAIR.gc3gc3.txt', header = TRUE)
r6 <- read.table('Rice.gc3gc3.txt', header = TRUE)

# Devide genes into Poor and Rich based on quantiles 10-90 
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

# Count % of monoexonic genes in poor and rich genes   ---------------------------
names(r1)[1] <- 'name'
names(r2)[1] <- 'name'
names(r3)[1] <- 'name'
names(r4)[1] <- 'name'

setwd('D:/larch/Annotations')

largff <- read.table('maker-res.gff', stringsAsFactors=FALSE, sep='\t', header=F, col.names = c('chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attribute'))
pitgff <- read.table('Pita.2_01.remake2.gff', stringsAsFactors=FALSE, sep='\t', header=F, col.names = c('chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attribute'))
pggff <- read.table('PG29v3-renamedID_1000nt.gff', stringsAsFactors=FALSE, sep='\t', header=F, col.names = c('chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attribute'))
abgff <- read.table('Pabies1.0.gff3', stringsAsFactors=FALSE, sep='\t', header=F, col.names = c('chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attribute'))

lar <- exons_merge(largff, r1)
pit <- exons_merge(pitgff, r4)
pab <- exons_merge(abgff, r2)
pg <- exons_merge(pggff, r3)

rm(largff, pitgff, pggff, abgff)

#  Summarise number of rows per gene name in merged table------------

pal <- c(poor="#3B4992B2", rich= "#EE0000B2")

fu <- as.data.frame(table(table(subset(lar, type=='poor')[,1])))
su <- as.data.frame(table(table(subset(lar, type=='rich')[,1])))
fu$type <- 'poor'; su$type <- 'rich' ;
nu <- rbind(fu, su) ; nu$species <- 'Larix'

      p1 <- ggplot()+
        geom_bar(data=nu, aes(x=Var1, y=Freq, fill=type), stat="identity", position = "dodge")+
        labs(title = 'L. sibirica', tag = 'a', x='Количество экзонов', y='Число генов')+
        theme_few()+
        theme(plot.title = element_text(hjust=0.5, face='italic'))+
        scale_fill_manual(values = pal)

fu <- as.data.frame(table(table(subset(pit, type=='poor')[,1])))
su <- as.data.frame(table(table(subset(pit, type=='rich')[,1])))
fu$type <- 'poor'; su$type <- 'rich' ;
nu <- rbind(fu, su) ; nu$species <- 'P.taeda'

      p2 <- ggplot()+
        geom_bar(data=nu, aes(x=Var1, y=Freq, fill=type), stat="identity", position = "dodge")+
        labs(title = 'P. taeda', tag = 'b', x='Количество экзонов', y='Число генов')+
        theme_few()+
        theme(plot.title = element_text(hjust=0.5, face='italic'))+
        scale_fill_manual(values = pal) 

fu <- as.data.frame(table(table(subset(pg, type=='poor')[,1])))
su <- as.data.frame(table(table(subset(pg, type=='rich')[,1])))
fu$type <- 'poor'; su$type <- 'rich' ;
nu <- rbind(fu, su) ; nu$species <- 'P.glauca'

      p3 <- ggplot()+
        geom_bar(data=nu, aes(x=Var1, y=Freq, fill=type), stat="identity", position = "dodge")+
        labs(title = 'P. glauca', tag = 'c', x='Количество экзонов', y='Число генов')+
        theme_few()+
        theme(plot.title = element_text(hjust=0.5, face='italic'))+
        scale_fill_manual(values = pal) 

fu <- as.data.frame(table(table(subset(pab, type=='poor')[,1]))) ; fu$type <- 'poor'
su <- as.data.frame(table(table(subset(pab, type=='rich')[,1]))) ; su$type <- 'rich'
nu <- rbind(fu, su) ; nu$species <- 'P.abies'

      p4 <- ggplot()+
        geom_bar(data=nu, aes(x=Var1, y=Freq, fill=type), stat="identity", position = "dodge")+
        labs(title = 'P. abies', tag = 'd', x='Количество экзонов', y='Число генов')+
        theme_few()+
        theme(plot.title = element_text(hjust=0.5, face='italic'))+
        scale_fill_manual(values = pal)
        





setwd('G:/Мой диск/Bioinform/Dissertation/Графики')

p <- (p1+p2)/(p3+p4) + plot_layout(guides = "collect") & theme(legend.position="bottom")

ggplot2::ggsave(filename = "23.Barplot.svg", 
                plot = p, 
                device = "svg", 
                width = 9,
                height = 6, 
                units = "in")
