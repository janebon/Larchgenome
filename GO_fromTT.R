
setwd('G:/Мой диск/Bioinform/Lar Annot/B2G/tsv files')

library(XNomial)
library(pwr)

custom <- c(Psme = 'darkseagreen',
            Potr = 'deepskyblue3', 
            Pita = 'darkseagreen2',
            Pila = 'darkseagreen3',
            Pgl = 'darkseagreen4',
            Pab = 'palegreen3',
            Lar = 'brown1',
            Bepe = 'dodgerblue4',
            Abal = 'darkseagreen1')


D = read.table('GO_table.txt', header = TRUE,  sep = "\t")
head(D)
domains = unique(D$domain)
BP = subset(D,D$domain == "Biological Process")
MF = subset(D,D$domain ==  "Molecular Function" )
CC = subset(D,D$domain == "Cellular Component")




TOTAL=read.table('Total_per_domain-1.txt',header=T,row.names = 1)
names(D)[1:12]=c('GO','domain','description','Psme','Potr','Pita','Pila','Pgl','Pab', 'Lar','Bepe','Abal')
TOTAL=TOTAL[c('Psme','Potr','Pita','Pila','Pgl','Pab', 'Lar','Bepe','Abal'),]
head(TOTAL)

# dir.create('GO_plots')
# setwd('GO_plots')

for(i in 1:dim(D)[1]){
  if(D$domain[i]=='Biological Process'){TOTAR=TOTAL[,1]; alpha=dim(BP)[1]}
  if(D$domain[i]=="Molecular Function"){TOTAR=TOTAL[,3]; alpha=dim(MF)[1]}
  if(D$domain[i]=="Cellular Component"){TOTAR=TOTAL[,2]; alpha=dim(CC)[1]}
  d=D[i,4:12]
  TD=sum(d)
  TT=sum(TOTAR)
  TOTAR=d/TT
  sum(TOTAR)
  # ch=chisq.test(d, p=TOTAR);ch
  # xmulti(obs=as.numeric(d),
  #        expr=TOTAR*TD,safety=1000000000,
  #       detail = 2)   
  # mult= xmonte(obs=as.numeric(d),  # Goodness-Of-Fit Test By Monte-Carlo
  #              expr=TOTAR*TD,
  #              detail = 2) 
  # 
  # p=pwr.chisq.test(w=ES.w1(TOTAR, as.numeric(d)/TD), N=TD, df=(length(d)-1), sig.level=0.01/alpha)   # power calculations for chi-squared tests
  # p$power
  
  # if(ch$p.value*alpha < 0.01 & p$power > 0.99){
    
    # write.table(paste(D$GO[i], D$description[i], ch$p.value, ch$p.value*alpha,
    #                   mult$asymptotic.p.value, p$power, sep='|'), 'STAT3.txt', append=T, row.names = FALSE, col.names = FALSE, quote=FALSE)
    # 
    # png(paste(gsub(":","_", D$GO[i]), 'png', sep='.'))
    # 
    # plot(as.numeric(d/TD), TOTAR, xlab='Observed', ylab='Expected', col=custom, pch=19, cex=2,
    #      main=paste(D$GO[i], D$description[i], 'p-value=', ch$p.value*alpha))
    # text(as.numeric(d/TD), TOTAR, labels=names(d), cex=1, font=2)
    # 
    # dev.off()
    # 
    
    if(D$domain[i]=='Biological Process'){dom='BP'}
    if(D$domain[i]=="Molecular Function"){dom='MF'}
    if(D$domain[i]=="Cellular Component"){dom='CC'}

    
    write.table(paste(D$GO[i], D$description[i], d/TD, TOTAR, names(d), dom, sep='\t'), 
                'STAT4.txt', append=T, row.names = FALSE, col.names = FALSE, quote=FALSE)
  # }
}



# Goodness-Of-Fit Test By Monte-Carlo
# power calculations for chi-squared tests

# parameters : ch$p.value * alpha < 0.01   &   p$power > 0.99




setwd('G:/Мой диск/Bioinform/Lar Annot/B2G/tsv files')
dat <- read.table('STAT4.txt', sep='\t')
names(dat) <- c('GO', 'description', 'Observed', 'Expected', 'species', 'domain')


library(ggplot2)

subset(dat, domain=='BP')
subset(dat, domain=='BP')

ggplot(dat, aes(x = Observed, y = Expected, color=species, shape=domain))+
  geom_point(size=4, alpha=0.5) +
  theme_bw() +
  scale_color_manual(values = custom)





