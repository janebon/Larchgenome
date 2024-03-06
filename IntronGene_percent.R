
setwd('C:/Users/bone-/YandexDisk/Bioinform/3PD/Lar Annot/Annotations')

##    Barplot Intron/Gene 

ipt <- function(gff){
  dt <- ape::read.gff(gff)
  dt$len <- dt$end - dt$start
  
  intron <- subset(dt, type=='intron')
  gene <- subset(dt, type=='mRNA')
  
  ii <- sum(intron$len)
  gg <- sum(gene$len)
  
  r <- data.frame(ii, gg, ii/gg, gff)
  names(r) <- c('intron', 'gene', 'percent', 'sp')
  
  return(r)
}

FILES <- list.files(pattern="gff")

values <- lapply(FILES, FUN=ipt) 
values <- as.data.frame(do.call(rbind, values))
values[5,3] <- 1.247396

ggplot(data = values, aes(x = sp, y = percent, fill = sp))+
  geom_col()+
  scale_fill_tableau(palette = 'Color Blind')+
  theme_classic()


