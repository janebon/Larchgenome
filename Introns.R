library(ggplot2)
library(ggthemes)
library(scales)
library(ape)
library(dplyr)
library(patchwork)

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



setwd('C:/Users/Пользователь/Google Диск/Bioinform/Lar Annot/Introns')


# Make intron length longer than 90q boxplot

intron.len.q <- function(file, qt){
  
  gff <- read.table(file, stringsAsFactors=FALSE, sep='\t', header=F, quote = '', comment.char = '#')
  names(gff) <- c('chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attribute')
  introns <- subset(gff, type == 'intron')
  introns$len <- introns$end - introns$start
  
  introns$name <- getAttributeField(introns, field='Parent')
  iu <- aggregate(introns$len, by=list(introns$name), FUN=sum)
  names(iu) <- c('name', 'len')
  
  fu <- subset(iu, len > quantile(iu$len, qt))
  fu$sp <- file
  return(fu)
  
}   # function to subset introns

FILES <- Sys.glob("*.*")
su <- lapply(FILES, FUN=intron.len.q, qt=0.9)  # subset introns longer than 90 quantile
nu <- as.data.frame(do.call(rbind, su))        # convert list of lists to data frame
nu <- subset(nu, sp!='O.sativa')
nu <- subset(nu, sp!='A.thaliana')

cols = c(
  L.sibirica ="firebrick1",
  
  A.alba = "olivedrab3",
  P.abies ="palegreen3",
  P.glauca ="seagreen4",
  P.taeda ="darkseagreen",
  P.lambertiana = "springgreen3",
  P.menziesii = "mediumaquamarine" ,
  P.tabuliformis = "mediumseagreen",
  
  F.sylvatica ="tan2",
  A.trichopoda = 'wheat', 
  P.trichocarpa ="peachpuff",
  V.vinifera = "lightgoldenrod",
  Z.mays = "navajowhite3",
  Q.robur = "orange")



level_order <- c('L.sibirica','A.alba', 'P.abies', 'P.glauca', 'P.taeda', 'P.lambertiana', 'P.menziesii', 'P.tabuliformis',
                 'A.trichopoda', 'F.sylvatica', 'P.trichocarpa', 'Q.robur', 'V.vinifera', 'Z.mays') 


format(data, decimal.mark=",")
options(OutDec= ",")

p <- ggplot(data = nu, aes(x=factor(sp, level = level_order), y=log10(len), col=sp))+
  geom_point(alpha = 0.4, size=1, position = "jitter")+
  geom_boxplot(width=0.6, outlier.shape = NA, color='black')+
  scale_y_continuous(breaks = pretty_breaks())+
  scale_color_manual(values = cols)+
  labs(x='Species', y=expression(bold(Log["10"]~"(intron length, bp)")))+
  theme_classic()+
  theme(legend.position = 'none',
        axis.text.x=element_text(size=16, face="italic", angle = 13, hjust = 0.5, vjust = 0.7), 
        axis.text.y=element_text(size=16))


setwd('C:/Users/Пользователь/Google Диск/Bioinform/Dissertation/Доп.мат/Plots')
ggplot2::ggsave(filename = "Introns.svg", 
                plot = p, 
                device = "svg", 
                width = 11,
                height =5, 
                dpi = 600,
                units = "in")
