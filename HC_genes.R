setwd('D:/')
setwd('/gpfs/HOME/ebondar/LarAnnot/blast_genes')
setwd('G:/Мой диск/Bioinform/Lar Annot/HC genes')


getAttributeField <- function(x, field, attrsep = ";") {
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



    # Load annotation and extract mRNAs
	
gff <- read.table("G:/Мой диск/Bioinform/Lar Annot/maker-res.wi.b2g.incmpl.renamed.gff", skip = 1, sep = '\t', quote=NULL, comment='', stringsAsFactors=FALSE)

gff <- read.table("/gpfs/HOME/ebondar/LarAnnot/maker-res.wi.b2g.incmpl.renamed.gff", skip = 1, sep = '\t', quote=NULL, comment='', stringsAsFactors=FALSE)

names(gff) <- c('chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attribute')

mrna <- subset(gff, type == 'mRNA')
mrna$id <- getAttributeField(mrna, 'ID')
mrna <- mrna[, c('chr', 'start', 'end', 'strand', 'id')]


    # Load RNA alignment to Genome
	
RNAgenome <- read.table('Lar_genome.count-2.txt', col.names = c('id', 'countG'))


    # Load RNA alignment to Genes
	
RNAgene <- read.table('Lar.genes.bedgraph', col.names = c('id', 'start', 'end', 'countg'), stringsAsFactors=FALSE)

RNAgene2 <- aggregate(RNAgene$countg, by=list(Category=RNAgene$id), FUN=max)
names(RNAgene2) <- c('id', 'countg')
RNAgene2$id <- as.character(RNAgene2$id)


    # Load BLAST results
	
blast <- read.table('maker-res.Embryophyta.clear.upd.psl', sep = '\t', quote=NULL, comment='')
names(blast) <- c('id', 'title', 'qlen', 'slen', 'length', 'evalue', 'pident')

blast <- blast[order(-blast$pident),] 
blast <- blast[!duplicated(blast[,'id']),]                         # remove duplicated genes
blast <- blast[, c('id', 'title', 'pident', 'qlen', 'length')]


    # Merge all evidences together

dt <- merge(mrna, blast, all = TRUE)  

dt <- merge(dt, RNAgenome)  ;  nrow(subset(dt, count > 10))

dt <- merge(dt, RNAgene2)


	# Filter
	
length(subset(dt, pident >= 70 & length/qlen*100 >= 50)$id)  # 13475 by blast hits
length(subset(dt, pident >= 60 & length/qlen*100 >= 30)$id)  # 23110 by blast hits
length(subset(dt, countG >= 5)$id)                           # 25002 by RNA alignment to Genome
length(subset(dt, countg >= 5)$id)                           # 30248 by RNA alignment to genes

length(subset(dt, pident >= 60 & length/qlen*100 >= 30 & countG >= 5)$id)  # 13716 by all
length(subset(dt, pident >= 70 & length/qlen*100 >= 30 & countG >= 5)$id)  # 8951
length(subset(dt, pident >= 70 & length/qlen*100 >= 30 & countG >= 3)$id)  # 9382
length(subset(dt, pident >= 30 & length/qlen*100 >= 50 & countG >= 3)$id)  # 21607

res <- subset(dt, pident >= 30 & length/qlen*100 >= 50 & countG >= 3)$id

write.table(res, 'HC_genes.list', quote=FALSE, row.names=FALSE, sep="\t")


# --------------
#  agat_sp_filter_feature_from_keep_list.pl --gff maker-res.wi.b2g.incmpl.renamed.gff --keep_list ./HC_genes/HC_genes.list  --output maker-res.wi.b2g.incmpl.renamed.HC.gff
 
  
gff <- read.table('D:/maker-res.wi.b2g.incmpl.renamed.HC.gff', skip = 1, sep = '\t', quote=NULL, comment='', stringsAsFactors=FALSE)
names(gff) <- c('chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attribute')

length(unique(gff$chr))

