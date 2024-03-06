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

library(data.table)

setwd('C:/Users/Пользователь/Google Диск/Bioinform/Lar Annot')

maker <- read.table('maker-res.wi.gff')                                         # load maker file
names(maker) <- c('chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attribute')

rna <- subset(maker, type=='mRNA')                                              # take RNAs from maker annotation
rna$id <- getAttributeField(rna, field='ID')                                    # parse attribute of maker RNAs
rna$parent <- getAttributeField(rna, field='Parent')
rna$aed <- getAttributeField(rna, field='_AED')
rna$eaed <- getAttributeField(rna, field='_eAED')
rna$attribute <- NULL


b2go <- fread('C:/Users/Пользователь/Google Диск/Bioinform/Lar Annot/B2G/blast2go_IntPro_blast_GO.gff3', skip=2)
names(b2go) <- c('id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attribute')
b2go$attribute <- gsub('ID=(.*)-mRNA-1;', '', b2go$attribute, perl=TRUE)        # drop extra ID info
b <- b2go[,c('id', 'attribute')]


res <- merge(rna, b, by='id')                                                   # merge maker & b2go 
res$attribute <- paste('ID=', res$id, ';Parent=', res$parent, ';_AED=', res$aed, ';_eAED=', res$eaed, ';', res$attribute, sep='')
res <- res[,-c(1,10,11,12)]

annot <- rbind(subset(maker, type!='mRNA'), res)                                # add merged RNAs to the rest of the maker annotation

annot$source <- ifelse(annot$type=='intron', 'GenomeTools', annot$source)       # add source for introns


write.table(annot, 'maker-res.wi.b2g.gff', quote=FALSE, row.names=FALSE, sep="\t")
